// Regression check for the CUDA backend of the photopeak PDFs.
//
// Builds a toy single-channel model and a four-channel simultaneous model from
// the library's PDFs, with 1M events per channel, and checks that:
//   1. the CPU and CUDA backends give the same NLL at many parameter points;
//   2. the CPU NLL at the start point still matches a recorded reference, which
//      catches value changes made by this library or by a ROOT bump, even when
//      both backends change together;
//   3. the CUDA backend, selected through BestAvailableBackend() as the
//      library's own fits do, beats the CPU backend by a minimum factor.
//
// It needs a GPU, so it cannot run inside the build sandbox. Run it with
// `nix run .#cuda-regression`; it exits non-zero if any check fails.

#include "RooFitUtils.hpp"

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

const Long64_t kEventsPerChannel = 1000000;
const Int_t kSimultaneousChannels = 4;
const Int_t kWarmupEvaluations = 5;
const Int_t kTimedEvaluations = 100;
const Int_t kAgreementPoints = 40;
const Double_t kMaxBackendRelativeDifference = 1e-10;
const Double_t kMaxReferenceRelativeDifference = 1e-9;

// Limits on the CUDA-to-CPU ratio of the median time per NLL evaluation. The
// simultaneous limit is looser because ROOT 6.40's CUDA path recomputes every
// channel on every evaluation, while the CPU path skips the unchanged ones.
// Tighten it once the pinned ROOT skips clean nodes on CUDA as well.
const Double_t kMaxSingleTimeRatio = 0.3;
const Double_t kMaxSimultaneousTimeRatio = 0.6;

// CPU NLL at the start point, recorded with ROOT 6.40.00. After a deliberate
// change to the model or the PDFs, replace these with the values the check
// prints.
const Double_t kReferenceSingleStartNll = -9884771.017917;
const Double_t kReferenceSimultaneousStartNll = -34051638.494087;

struct ToyChannel {
  TString name;
  Double_t center;
  RooAddPdf *pdf;
  RooDataSet *data;
};

struct ToyModel {
  TString label;
  RooAbsPdf *pdf;
  RooAbsData *data;
  RooArgList floating;
  Bool_t simultaneous;
};

ToyChannel BuildChannel(Int_t index, RooRealVar &x, RooRealVar &lowTau,
                        RooRealVar &highTau, RooArgList &floating) {
  TString name = TString::Format("ch%d", index);
  Double_t center = 40.0 + 40.0 * index;
  x.setRange("fitrange_" + name, center - 20.0, center + 20.0);

  RooRealVar *mu = new RooRealVar("mu_" + name, "", center + 0.3, center - 5.0,
                                  center + 5.0);
  RooRealVar *sigma = new RooRealVar("sigma_" + name, "", 1.8, 0.5, 5.0);
  RooRealVar *linSlope = new RooRealVar("linSlope_" + name, "", 0.05, 0.0, 0.5);
  RooRealVar *bkgSlope =
      new RooRealVar("bkgSlope_" + name, "", 0.0, -0.004, 0.004);

  RooArgList components;
  components.add(*new RooGaussian("gauss_" + name, "", x, *mu, *sigma));
  components.add(
      *new RooLowExpTail("lowExp_" + name, "", x, *mu, *sigma, lowTau));
  components.add(*new RooStepShelf("shelf_" + name, "", x, *mu, *sigma));
  components.add(
      *new RooLowLinTail("lowLin_" + name, "", x, *mu, *sigma, *linSlope));
  components.add(
      *new RooHighExpTail("highExp_" + name, "", x, *mu, *sigma, highTau));
  components.add(
      *new RooPolynomial("bkg_" + name, "", x, RooArgList(*bkgSlope)));

  const Double_t fractions[6] = {0.60, 0.15, 0.10, 0.05, 0.05, 0.05};
  const char *yieldPrefixes[6] = {"nGauss_",  "nLowExp_",  "nShelf_",
                                  "nLowLin_", "nHighExp_", "nBkg_"};
  RooArgList yields;
  for (Int_t k = 0; k < 6; ++k) {
    RooRealVar *yield = new RooRealVar(yieldPrefixes[k] + name, "",
                                       fractions[k] * kEventsPerChannel * 1.02,
                                       0.0, 3.0 * kEventsPerChannel);
    yields.add(*yield);
    floating.add(*yield);
  }
  floating.add(RooArgList(*mu, *sigma, *linSlope, *bkgSlope));

  ToyChannel channel;
  channel.name = name;
  channel.center = center;
  channel.pdf = new RooAddPdf("model_" + name, "", components, yields);

  // Shapes roughly like the model's components, so the fit is well posed. The
  // seed fixes the data, which the reference NLL values depend on.
  TRandom3 rng(1000 + index);
  channel.data = new RooDataSet("data_" + name, "", RooArgSet(x));
  Long64_t accepted = 0;
  while (accepted < kEventsPerChannel) {
    Double_t r = rng.Uniform();
    Double_t value;
    if (r < 0.60)
      value = center + 2.0 * rng.Gaus();
    else if (r < 0.75)
      value = center - rng.Exp(6.0) + 2.0 * rng.Gaus();
    else if (r < 0.85)
      value = rng.Uniform(center - 20.0, center);
    else if (r < 0.90)
      value = center - rng.Uniform(0.0, 15.0);
    else if (r < 0.95)
      value = center + rng.Exp(4.0) + 2.0 * rng.Gaus();
    else
      value = rng.Uniform(center - 20.0, center + 20.0);
    if (value < center - 20.0 || value > center + 20.0)
      continue;
    x.setVal(value);
    channel.data->add(RooArgSet(x));
    ++accepted;
  }
  return channel;
}

ToyModel BuildModel(Int_t nChannels) {
  ToyModel model;
  model.simultaneous = nChannels > 1;
  model.label = model.simultaneous
                    ? TString::Format("simultaneous, %d x %lld events",
                                      nChannels, kEventsPerChannel)
                    : TString::Format("single, %lld events", kEventsPerChannel);

  RooRealVar *x = new RooRealVar("x", "x", 0.0, 40.0 + 40.0 * nChannels);
  RooRealVar *lowTau = new RooRealVar("lowTau", "", 3.2, 0.5, 10.0);
  RooRealVar *highTau = new RooRealVar("highTau", "", 2.1, 0.5, 10.0);
  model.floating.add(RooArgList(*lowTau, *highTau));

  std::vector<ToyChannel> channels;
  for (Int_t i = 0; i < nChannels; ++i)
    channels.push_back(BuildChannel(i, *x, *lowTau, *highTau, model.floating));
  Double_t low = channels.front().center - 20.0;
  Double_t high = channels.back().center + 20.0;
  x->setRange("fitrange", low, high);

  if (!model.simultaneous) {
    model.pdf = channels[0].pdf;
    model.data = channels[0].data;
    return model;
  }
  RooCategory *category = new RooCategory("channel", "");
  RooSimultaneous *sim = new RooSimultaneous("sim", "", *category);
  std::map<std::string, RooDataSet *> perChannel;
  for (size_t i = 0; i < channels.size(); ++i) {
    // As RooFitUtils::FitSimultaneous does before a SplitRange fit.
    channels[i].pdf->fixCoefRange("fitrange_" + channels[i].name);
    category->defineType(channels[i].name.Data());
    sim->addPdf(*channels[i].pdf, channels[i].name.Data());
    perChannel[channels[i].name.Data()] = channels[i].data;
  }
  model.pdf = sim;
  model.data =
      new RooDataSet("combined", "", RooArgSet(*x), RooFit::Index(*category),
                     RooFit::Import(perChannel));
  return model;
}

RooAbsReal *CreateNll(ToyModel &model, const RooCmdArg &backend) {
  if (model.simultaneous)
    return model.pdf->createNLL(*model.data, backend, RooFit::Extended(kTRUE),
                                RooFit::Range("fitrange"),
                                RooFit::SplitRange(kTRUE));
  return model.pdf->createNLL(*model.data, backend, RooFit::Extended(kTRUE),
                              RooFit::Range("fitrange"));
}

Double_t MaxRelativeDifference(RooAbsReal &cpuNll, RooAbsReal &cudaNll,
                               RooArgList &floating) {
  std::vector<Double_t> start;
  for (Int_t i = 0; i < floating.getSize(); ++i)
    start.push_back(static_cast<RooRealVar *>(floating.at(i))->getVal());

  TRandom3 rng(4242);
  Double_t worst = 0.0;
  for (Int_t point = 0; point < kAgreementPoints; ++point) {
    for (Int_t i = 0; i < floating.getSize(); ++i)
      static_cast<RooRealVar *>(floating.at(i))
          ->setVal(start[i] * (1.0 + 0.01 * rng.Uniform(-1.0, 1.0)));
    Double_t cpu = cpuNll.getVal();
    Double_t cuda = cudaNll.getVal();
    if (!std::isfinite(cpu) || !std::isfinite(cuda)) {
      worst = std::numeric_limits<Double_t>::infinity();
      break;
    }
    worst =
        std::max(worst, std::abs(cpu - cuda) / std::max(1.0, std::abs(cpu)));
  }

  for (Int_t i = 0; i < floating.getSize(); ++i)
    static_cast<RooRealVar *>(floating.at(i))->setVal(start[i]);
  return worst;
}

// Each evaluation moves one parameter, as Migrad's numerical gradient does, so
// a backend that skips unchanged nodes gets the benefit it would get in a fit.
Double_t MedianEvaluationMillis(RooAbsReal &nll, RooArgList &floating) {
  for (Int_t i = 0; i < kWarmupEvaluations; ++i)
    nll.getVal();

  std::vector<Double_t> millis;
  for (Int_t i = 0; i < kTimedEvaluations; ++i) {
    RooRealVar *parameter =
        static_cast<RooRealVar *>(floating.at(i % floating.getSize()));
    Double_t original = parameter->getVal();
    parameter->setVal(original + 1e-4 * (std::abs(original) + 1e-3));
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    nll.getVal();
    std::chrono::steady_clock::time_point stop =
        std::chrono::steady_clock::now();
    millis.push_back(
        std::chrono::duration<Double_t, std::milli>(stop - start).count());
    parameter->setVal(original);
  }
  std::sort(millis.begin(), millis.end());
  return millis[millis.size() / 2];
}

Bool_t Report(const TString &check, Bool_t passed, const TString &detail) {
  std::cout << (passed ? "  PASS  " : "  FAIL  ") << check << ": " << detail
            << std::endl;
  return passed;
}

Bool_t CheckModel(Int_t nChannels, Double_t referenceStartNll,
                  Double_t maxTimeRatio) {
  ToyModel model = BuildModel(nChannels);
  std::cout << model.label << " (" << model.floating.getSize()
            << " floating parameters)" << std::endl;

  RooAbsReal *cpuNll = CreateNll(model, RooFit::EvalBackend::Cpu());
  RooAbsReal *cudaNll = CreateNll(model, BestAvailableBackend());

  Double_t startNll = cpuNll->getVal();
  Double_t referenceDifference = std::abs(startNll - referenceStartNll) /
                                 std::max(1.0, std::abs(referenceStartNll));
  Double_t backendDifference =
      MaxRelativeDifference(*cpuNll, *cudaNll, model.floating);
  Double_t cpuMillis = MedianEvaluationMillis(*cpuNll, model.floating);
  Double_t cudaMillis = MedianEvaluationMillis(*cudaNll, model.floating);
  Double_t timeRatio = cudaMillis / cpuMillis;

  Bool_t passed = kTRUE;
  passed &= Report(
      "CPU and CUDA agree", backendDifference <= kMaxBackendRelativeDifference,
      TString::Format(
          "max relative NLL difference %.2e over %d points (limit %.0e)",
          backendDifference, kAgreementPoints, kMaxBackendRelativeDifference));
  passed &= Report(
      "start NLL matches reference",
      referenceDifference <= kMaxReferenceRelativeDifference,
      TString::Format(
          "%.6f vs reference %.6f (relative difference %.2e, limit %.0e)",
          startNll, referenceStartNll, referenceDifference,
          kMaxReferenceRelativeDifference));
  passed &= Report(
      "CUDA is fast enough", timeRatio <= maxTimeRatio,
      TString::Format("median per evaluation: CUDA %.3f ms, CPU %.3f ms, "
                      "ratio %.3f (limit %.2f)",
                      cudaMillis, cpuMillis, timeRatio, maxTimeRatio));
  delete cpuNll;
  delete cudaNll;
  return passed;
}

} // namespace

int main() {
#if !defined(AU_ROOFIT_BACKEND_CUDA) || !AU_ROOFIT_BACKEND_CUDA
  std::cout << "cuda-regression must be built with -DAU_ROOFIT_BACKEND_CUDA=1"
            << std::endl;
  return 2;
#else
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  std::cout << "CUDA backend regression check, ROOT " << gROOT->GetVersion()
            << std::endl;

  Bool_t passed = kTRUE;
  passed &= CheckModel(1, kReferenceSingleStartNll, kMaxSingleTimeRatio);
  passed &= CheckModel(kSimultaneousChannels, kReferenceSimultaneousStartNll,
                       kMaxSimultaneousTimeRatio);

  std::cout << (passed ? "All checks passed." : "Some checks FAILED.")
            << std::endl;
  return passed ? 0 : 1;
#endif
}
