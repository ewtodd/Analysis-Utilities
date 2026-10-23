// Regression check for the CUDA backend and for the likelihood RooFit builds.
//
// Builds toy models from the library's PDFs, with 1M events per channel, in the
// three configurations the library fits: a single channel, four channels with
// SplitRange (RooFitUtils::FitSimultaneous) and four channels over one shared
// range (the interactive simultaneous editor). As in the library, each
// exponential tail's yield is the Gaussian yield times a full-range ratio times
// a RooExpTailWindowFactor, a scalar node evaluated on the CPU in both
// backends. For each configuration it checks that:
//   1. the CPU and CUDA backends give the same NLL at many parameter points;
//   2. RooFit's NLL equals the likelihood the model defines, computed directly
//      from the components' values and integrals, up to a constant, at several
//      parameter points (released ROOT up to 6.40.04 fails this for the
//      simultaneous configurations, root-project/root#23444);
//   3. the CPU NLL at the start point still matches a recorded reference, which
//      also catches changes to the PDFs themselves;
//   4. the CUDA backend, selected through BestAvailableBackend() as the
//      library's own fits do, beats the CPU backend by a minimum factor (single
//      and SplitRange configurations);
//   5. Migrad fits with likelihood offsetting, which every fit in the library
//      uses, converge and reach the same parameters on CPU and CUDA (single
//      and SplitRange configurations). RooFit applies the offset only inside
//      the minimiser, so only a fit exercises it. These fits use a simpler
//      model generated from itself, with 100k events per channel, so that each
//      has one clear minimum and the CPU fits stay short.
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
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGlobalFunc.h>
#include <RooMinimizer.h>
#include <RooMsgService.h>
#include <RooPolynomial.h>
#include <RooRandom.h>
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
const Long64_t kFitEventsPerChannel = 100000;
const Int_t kSimultaneousChannels = 4;
const Int_t kWarmupEvaluations = 5;
const Int_t kTimedEvaluations = 100;
const Int_t kAgreementPoints = 40;
const Int_t kDirectPoints = 3;
const Double_t kMaxBackendRelativeDifference = 1e-10;
const Double_t kMaxDirectRelativeMismatch = 1e-9;
const Double_t kMaxReferenceRelativeDifference = 1e-9;
// Largest CPU/CUDA difference of a fitted parameter, in units of its error.
// Each fit stops somewhere within an EDM below 1e-3 of the minimum, so a few
// hundredths of an error are expected; a broken offset path is far outside.
const Double_t kMaxFitSigmaDifference = 0.05;
const Double_t kMaxConvergedEdm = 1e-3;

// Limit on the CUDA-to-CPU ratio of the median time per NLL evaluation.
const Double_t kMaxTimeRatio = 0.15;

// CPU NLL at the start point, recorded with ROOT master e0e5723692 after the
// direct-likelihood check passed. After a deliberate change to the model or
// the PDFs, replace these with the values the check prints. The shared-range
// value was re-recorded when the exponential tails moved to full-range ratios:
// the start values of those ratios reproduce the old tail yields in each
// channel's own window, so the other two did not change, but over the shared
// range a fixed ratio puts more of each tail inside the window.
const Double_t kReferenceSingleStartNll = -9903979.467281;
const Double_t kReferenceSplitStartNll = -34073745.226293;
const Double_t kReferenceSharedStartNll = -33521051.546687;

struct ToyChannel {
  TString name;
  Double_t center;
  RooAddPdf *pdf;
  RooDataSet *data;
  TString yieldRange;    ///< The range this channel's yields count events in.
  RooRealVar *windowLow; ///< yieldRange's edges, for the window factors.
  RooRealVar *windowHigh;
};

struct ToyModel {
  RooRealVar *x;
  std::vector<ToyChannel> channels;
  RooAbsPdf *pdf;
  RooAbsData *data;
  RooArgList floating;
  Bool_t simultaneous;
};

// As RooFitUtils builds it: an exponential tail's yield is nGauss * R times a
// RooExpTailWindowFactor, with R the full-range tail/Gaussian yield ratio and
// the factor taking it into the channel's window. R starts at the value that
// gives the tail startYield events, so the starting model is the one a plain
// yield of startYield defines.
RooAbsReal *ExpTailYield(const TString &prefix, ToyChannel &channel,
                         RooRealVar &nGauss, RooRealVar &mu, RooRealVar &sigma,
                         RooRealVar &tau, Bool_t highSide, Double_t startYield,
                         RooArgList &floating) {
  RooExpTailWindowFactor *factor = new RooExpTailWindowFactor(
      "factor" + prefix + channel.name, "", mu, sigma, tau, *channel.windowLow,
      *channel.windowHigh, highSide);
  RooRealVar *ratio =
      new RooRealVar("r" + prefix + channel.name, "",
                     startYield / (nGauss.getVal() * factor->getVal()), 0.0,
                     RooFitUtils::kExpTailAmplitudeMax);
  floating.add(*ratio);
  return new RooFormulaVar("n" + prefix + channel.name, "@0*@1*@2",
                           RooArgList(nGauss, *ratio, *factor));
}

ToyChannel BuildChannel(Int_t index, Long64_t events, RooRealVar &x,
                        RooRealVar &lowTau, RooRealVar &highTau,
                        RooArgList &floating) {
  TString name = TString::Format("ch%d", index);
  Double_t center = 40.0 + 40.0 * index;
  x.setRange("fitrange_" + name, center - 20.0, center + 20.0);

  ToyChannel channel;
  channel.name = name;
  channel.center = center;
  channel.windowLow = new RooRealVar("windowLow_" + name, "", center - 20.0);
  channel.windowHigh = new RooRealVar("windowHigh_" + name, "", center + 20.0);

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
  RooRealVar *nGauss =
      new RooRealVar(yieldPrefixes[0] + name, "", fractions[0] * events * 1.02,
                     0.0, 3.0 * events);
  floating.add(*nGauss);
  RooArgList yields;
  yields.add(*nGauss);
  for (Int_t k = 1; k < 6; ++k) {
    Double_t startYield = fractions[k] * events * 1.02;
    if (k == 1) {
      yields.add(*ExpTailYield("LowExp_", channel, *nGauss, *mu, *sigma, lowTau,
                               kFALSE, startYield, floating));
    } else if (k == 4) {
      yields.add(*ExpTailYield("HighExp_", channel, *nGauss, *mu, *sigma,
                               highTau, kTRUE, startYield, floating));
    } else {
      RooRealVar *yield = new RooRealVar(yieldPrefixes[k] + name, "",
                                         startYield, 0.0, 3.0 * events);
      yields.add(*yield);
      floating.add(*yield);
    }
  }
  floating.add(RooArgList(*mu, *sigma, *linSlope, *bkgSlope));

  channel.pdf = new RooAddPdf("model_" + name, "", components, yields);

  // Shapes roughly like the model's components, so the fit is well posed. The
  // seed fixes the data, which the reference NLL values depend on.
  TRandom3 rng(1000 + index);
  channel.data = new RooDataSet("data_" + name, "", RooArgSet(x));
  Long64_t accepted = 0;
  while (accepted < events) {
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

// As RooFitUtils does: the default range of x is the fit range (one channel)
// or the union of the channel windows (several), under the name "fitrange".
void FinishModel(ToyModel &model) {
  Double_t low = model.channels.front().center - 20.0;
  Double_t high = model.channels.back().center + 20.0;
  model.x->setRange(low, high);
  model.x->setRange("fitrange", low, high);

  if (!model.simultaneous) {
    model.channels[0].yieldRange = "fitrange";
    model.pdf = model.channels[0].pdf;
    model.data = model.channels[0].data;
    return;
  }
  RooCategory *category = new RooCategory("channel", "");
  RooSimultaneous *sim = new RooSimultaneous("sim", "", *category);
  std::map<std::string, RooDataSet *> perChannel;
  for (size_t i = 0; i < model.channels.size(); ++i) {
    category->defineType(model.channels[i].name.Data());
    sim->addPdf(*model.channels[i].pdf, model.channels[i].name.Data());
    perChannel[model.channels[i].name.Data()] = model.channels[i].data;
  }
  model.pdf = sim;
  model.data =
      new RooDataSet("combined", "", RooArgSet(*model.x),
                     RooFit::Index(*category), RooFit::Import(perChannel));
}

ToyModel BuildModel(Int_t nChannels, Long64_t events) {
  ToyModel model;
  model.simultaneous = nChannels > 1;
  model.x = new RooRealVar("x", "x", 0.0, 40.0 + 40.0 * nChannels);
  RooRealVar *lowTau = new RooRealVar("lowTau", "", 3.2, 0.5, 10.0);
  RooRealVar *highTau = new RooRealVar("highTau", "", 2.1, 0.5, 10.0);
  model.floating.add(RooArgList(*lowTau, *highTau));

  for (Int_t i = 0; i < nChannels; ++i)
    model.channels.push_back(
        BuildChannel(i, events, *model.x, *lowTau, *highTau, model.floating));
  FinishModel(model);
  return model;
}

// The model for the fit checks: per channel a Gaussian, both exponential tails
// and a linear background, with data generated from the model itself and every
// parameter well inside its limits, so that a fit has one clear minimum. The
// fit then starts from slightly displaced values.
ToyModel BuildFitModel(Int_t nChannels, Long64_t events) {
  RooRandom::randomGenerator()->SetSeed(2024);
  ToyModel model;
  model.simultaneous = nChannels > 1;
  model.x = new RooRealVar("x", "x", 0.0, 40.0 + 40.0 * nChannels);
  RooRealVar *lowTau = new RooRealVar("lowTau", "", 3.0, 0.5, 10.0);
  RooRealVar *highTau = new RooRealVar("highTau", "", 2.0, 0.5, 10.0);
  model.floating.add(RooArgList(*lowTau, *highTau));

  const Double_t fractions[4] = {0.75, 0.10, 0.05, 0.10};
  const char *yieldPrefixes[4] = {"nGauss_", "nLowExp_", "nHighExp_", "nBkg_"};
  for (Int_t i = 0; i < nChannels; ++i) {
    ToyChannel channel;
    channel.name = TString::Format("ch%d", i);
    channel.center = 40.0 + 40.0 * i;
    Double_t low = channel.center - 20.0;
    Double_t high = channel.center + 20.0;
    model.x->setRange("fitrange_" + channel.name, low, high);
    channel.windowLow = new RooRealVar("windowLow_" + channel.name, "", low);
    channel.windowHigh = new RooRealVar("windowHigh_" + channel.name, "", high);

    RooRealVar *mu = new RooRealVar("mu_" + channel.name, "", channel.center,
                                    channel.center - 5.0, channel.center + 5.0);
    RooRealVar *sigma =
        new RooRealVar("sigma_" + channel.name, "", 1.5, 0.5, 5.0);
    RooRealVar *bkgSlope =
        new RooRealVar("bkgSlope_" + channel.name, "", 0.002, -0.004, 0.004);
    RooArgList components;
    components.add(
        *new RooGaussian("gauss_" + channel.name, "", *model.x, *mu, *sigma));
    components.add(*new RooLowExpTail("lowExp_" + channel.name, "", *model.x,
                                      *mu, *sigma, *lowTau));
    components.add(*new RooHighExpTail("highExp_" + channel.name, "", *model.x,
                                       *mu, *sigma, *highTau));
    components.add(*new RooPolynomial("bkg_" + channel.name, "", *model.x,
                                      RooArgList(*bkgSlope)));
    RooRealVar *nGauss =
        new RooRealVar(yieldPrefixes[0] + channel.name, "",
                       fractions[0] * events, 0.0, 3.0 * events);
    RooRealVar *nBkg = new RooRealVar(yieldPrefixes[3] + channel.name, "",
                                      fractions[3] * events, 0.0, 3.0 * events);
    model.floating.add(*nGauss);
    RooArgList yields;
    yields.add(*nGauss);
    yields.add(*ExpTailYield("LowExp_", channel, *nGauss, *mu, *sigma, *lowTau,
                             kFALSE, fractions[1] * events, model.floating));
    yields.add(*ExpTailYield("HighExp_", channel, *nGauss, *mu, *sigma,
                             *highTau, kTRUE, fractions[2] * events,
                             model.floating));
    yields.add(*nBkg);
    model.floating.add(*nBkg);
    model.floating.add(RooArgList(*mu, *sigma, *bkgSlope));
    channel.pdf =
        new RooAddPdf("model_" + channel.name, "", components, yields);

    // Generated over the channel's window, which is x's default range here.
    model.x->setRange(low, high);
    channel.data = channel.pdf->generate(RooArgSet(*model.x), events);
    model.channels.push_back(channel);
  }
  FinishModel(model);

  for (Int_t i = 0; i < model.floating.getSize(); ++i) {
    RooRealVar *parameter = static_cast<RooRealVar *>(model.floating.at(i));
    TString name = parameter->GetName();
    if (name.BeginsWith("n"))
      parameter->setVal(parameter->getVal() * 1.03);
    else if (name.BeginsWith("mu_"))
      parameter->setVal(parameter->getVal() + 0.15);
    else if (name.BeginsWith("sigma_"))
      parameter->setVal(parameter->getVal() * 1.05);
    else if (name.BeginsWith("bkgSlope_"))
      parameter->setVal(0.001);
    else if (name == "lowTau")
      parameter->setVal(3.3);
    else if (name == "highTau")
      parameter->setVal(1.8);
  }
  return model;
}

// As RooFitUtils::FixChannelCoefRanges(): each channel's own window for a
// SplitRange fit, the shared "fitrange" otherwise, with the window factors
// moved to the same range.
void SetYieldRanges(ToyModel &model, Bool_t splitRange) {
  for (size_t i = 0; i < model.channels.size(); ++i) {
    ToyChannel &channel = model.channels[i];
    channel.yieldRange =
        splitRange ? TString("fitrange_" + channel.name) : TString("fitrange");
    channel.pdf->fixCoefRange(channel.yieldRange);
    channel.windowLow->setVal(model.x->getMin(channel.yieldRange));
    channel.windowHigh->setVal(model.x->getMax(channel.yieldRange));
  }
}

RooAbsReal *CreateNll(ToyModel &model, Bool_t splitRange,
                      const RooCmdArg &backend, Bool_t offset) {
  RooCmdArg offsetArg = RooFit::Offset(offset ? "initial" : "off");
  if (splitRange)
    return model.pdf->createNLL(*model.data, backend, RooFit::Extended(kTRUE),
                                RooFit::Range("fitrange"),
                                RooFit::SplitRange(kTRUE), offsetArg);
  return model.pdf->createNLL(*model.data, backend, RooFit::Extended(kTRUE),
                              RooFit::Range("fitrange"), offsetArg);
}

// The extended NLL the model defines, computed without RooFit's likelihood
// machinery: for each channel, the sum of its yields y_k minus, over its
// events, log(sum_k y_k f_k(x) / I_k), where f_k is component k's unnormalised
// value and I_k its integral over the range the channel's yields refer to.
Double_t DirectNll(ToyModel &model) {
  Double_t nll = 0.0;
  for (size_t c = 0; c < model.channels.size(); ++c) {
    ToyChannel &channel = model.channels[c];
    RooArgList components(channel.pdf->pdfList());
    RooArgList yields(channel.pdf->coefList());
    std::vector<Double_t> integral;
    std::vector<Double_t> yield;
    for (Int_t k = 0; k < components.getSize(); ++k) {
      RooAbsReal *integralOverRange =
          static_cast<RooAbsPdf *>(components.at(k))
              ->createIntegral(RooArgSet(*model.x),
                               RooFit::Range(channel.yieldRange));
      integral.push_back(integralOverRange->getVal());
      delete integralOverRange;
      yield.push_back(static_cast<RooAbsReal *>(yields.at(k))->getVal());
      nll += yield.back();
    }
    for (Int_t i = 0; i < channel.data->numEntries(); ++i) {
      model.x->setVal(channel.data->get(i)->getRealValue("x"));
      Double_t density = 0.0;
      for (Int_t k = 0; k < components.getSize(); ++k)
        density += yield[k] *
                   static_cast<RooAbsPdf *>(components.at(k))->getVal() /
                   integral[k];
      nll -= std::log(density);
    }
  }
  return nll;
}

// Point 0 is the start; the others move yields, shapes and positions.
void MoveToPoint(RooArgList &floating, const std::vector<Double_t> &start,
                 Int_t point) {
  for (Int_t i = 0; i < floating.getSize(); ++i) {
    RooRealVar *parameter = static_cast<RooRealVar *>(floating.at(i));
    TString name = parameter->GetName();
    Double_t value = start[i];
    if (point == 1) {
      if (name.BeginsWith("n"))
        value *= 1.05;
      else if (name.BeginsWith("sigma_"))
        value *= 1.1;
      else if (name == "lowTau")
        value += 0.3;
    } else if (point == 2) {
      if (name.BeginsWith("n"))
        value *= 0.97;
      else if (name.BeginsWith("mu_"))
        value += 0.2;
      else if (name == "highTau")
        value -= 0.2;
      else if (name.BeginsWith("linSlope_"))
        value *= 1.2;
      else if (name.BeginsWith("rLowExp_") || name.BeginsWith("rHighExp_"))
        value *= 1.1;
    }
    parameter->setVal(value);
  }
}

std::vector<Double_t> Values(RooArgList &floating) {
  std::vector<Double_t> values;
  for (Int_t i = 0; i < floating.getSize(); ++i)
    values.push_back(static_cast<RooRealVar *>(floating.at(i))->getVal());
  return values;
}

// RooFit adds parameter-independent terms (one for the category of a
// simultaneous model, for example), so compare how (RooFit - direct) changes
// between parameter points rather than its value.
Double_t MaxDirectMismatch(ToyModel &model, RooAbsReal &nll) {
  std::vector<Double_t> start = Values(model.floating);
  Double_t reference = 0.0;
  Double_t worst = 0.0;
  for (Int_t point = 0; point < kDirectPoints; ++point) {
    MoveToPoint(model.floating, start, point);
    Double_t difference = nll.getVal() - DirectNll(model);
    if (!std::isfinite(difference)) {
      worst = std::numeric_limits<Double_t>::infinity();
      break;
    }
    if (point == 0)
      reference = difference;
    worst = std::max(worst, std::abs(difference - reference));
  }
  MoveToPoint(model.floating, start, 0);
  return worst;
}

// The largest absolute difference; the caller scales it by the size of the NLL.
Double_t MaxBackendDifference(RooAbsReal &cpuNll, RooAbsReal &cudaNll,
                              RooArgList &floating) {
  std::vector<Double_t> start = Values(floating);
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
    worst = std::max(worst, std::abs(cpu - cuda));
  }
  MoveToPoint(floating, start, 0);
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

Bool_t CheckConfiguration(ToyModel &model, const TString &label,
                          Bool_t splitRange, Double_t referenceStartNll,
                          Bool_t checkSpeed) {
  std::cout << label << " (" << model.floating.getSize()
            << " floating parameters)" << std::endl;
  RooAbsReal *cpuNll =
      CreateNll(model, splitRange, RooFit::EvalBackend::Cpu(), kFALSE);
  RooAbsReal *cudaNll =
      CreateNll(model, splitRange, BestAvailableBackend(), kFALSE);

  Double_t startNll = cpuNll->getVal();
  Double_t backendDifference =
      MaxBackendDifference(*cpuNll, *cudaNll, model.floating) /
      std::max(1.0, std::abs(startNll));
  Double_t directMismatch = MaxDirectMismatch(model, *cpuNll);
  Double_t directLimit = kMaxDirectRelativeMismatch * std::abs(startNll);
  Double_t referenceDifference = std::abs(startNll - referenceStartNll) /
                                 std::max(1.0, std::abs(referenceStartNll));

  Bool_t passed = kTRUE;
  passed &= Report(
      "CPU and CUDA agree", backendDifference <= kMaxBackendRelativeDifference,
      TString::Format(
          "max NLL difference %.2e of the NLL over %d points (limit %.0e)",
          backendDifference, kAgreementPoints, kMaxBackendRelativeDifference));
  passed &= Report(
      "NLL matches the direct likelihood", directMismatch <= directLimit,
      TString::Format("(RooFit - direct) varies by %.3g over %d parameter "
                      "points (limit %.3g)",
                      directMismatch, kDirectPoints, directLimit));
  passed &= Report(
      "start NLL matches reference",
      referenceDifference <= kMaxReferenceRelativeDifference,
      TString::Format(
          "%.6f vs reference %.6f (relative difference %.2e, limit %.0e)",
          startNll, referenceStartNll, referenceDifference,
          kMaxReferenceRelativeDifference));
  if (checkSpeed) {
    Double_t cpuMillis = MedianEvaluationMillis(*cpuNll, model.floating);
    Double_t cudaMillis = MedianEvaluationMillis(*cudaNll, model.floating);
    Double_t timeRatio = cudaMillis / cpuMillis;
    passed &= Report(
        "CUDA is fast enough", timeRatio <= kMaxTimeRatio,
        TString::Format("median per evaluation: CUDA %.3f ms, CPU %.3f ms, "
                        "ratio %.3f (limit %.2f)",
                        cudaMillis, cpuMillis, timeRatio, kMaxTimeRatio));
  }
  delete cpuNll;
  delete cudaNll;
  return passed;
}

struct OffsetFit {
  Int_t status;
  Double_t edm;
  Double_t startNll; ///< Absolute NLL at the starting point.
  Double_t minNll;   ///< As Minuit saw it, so relative to the start if offset.
  std::vector<Double_t> values;
  std::vector<Double_t> errors;
};

// Configured as RooFitUtils' fits are: Minuit2 Migrad with
// RooFit::Offset("initial"). The parameters are put back afterwards.
OffsetFit FitWithOffset(ToyModel &model, Bool_t splitRange,
                        const RooCmdArg &backend) {
  std::vector<Double_t> start = Values(model.floating);
  RooAbsReal *nll = CreateNll(model, splitRange, backend, kTRUE);
  RooMinimizer minimizer(*nll);
  minimizer.setPrintLevel(-1);
  minimizer.setStrategy(1);
  minimizer.setOffsetting(kTRUE);

  OffsetFit fit;
  fit.startNll = nll->getVal();
  fit.status = minimizer.minimize("Minuit2", "migrad");
  RooFitResult *result = minimizer.save();
  fit.edm = result->edm();
  fit.minNll = result->minNll();
  delete result;
  fit.values = Values(model.floating);
  for (Int_t i = 0; i < model.floating.getSize(); ++i)
    fit.errors.push_back(
        static_cast<RooRealVar *>(model.floating.at(i))->getError());

  MoveToPoint(model.floating, start, 0);
  delete nll;
  return fit;
}

Bool_t CheckOffsetFits(ToyModel &model, const TString &label,
                       Bool_t splitRange) {
  OffsetFit cpu = FitWithOffset(model, splitRange, RooFit::EvalBackend::Cpu());
  OffsetFit cuda = FitWithOffset(model, splitRange, BestAvailableBackend());

  Double_t worst = 0.0;
  for (Int_t i = 0; i < model.floating.getSize(); ++i)
    worst = std::max(worst,
                     std::abs(cpu.values[i] - cuda.values[i]) / cpu.errors[i]);
  Bool_t converged = cpu.status == 0 && cuda.status == 0 &&
                     cpu.edm < kMaxConvergedEdm && cuda.edm < kMaxConvergedEdm;
  // With offsetting, Minuit's minimum is the improvement over the start, a
  // tiny fraction of the absolute NLL. Without it the comparison would not
  // test the offset path at all.
  Bool_t offsetUsed = std::abs(cpu.minNll) < 1e-3 * std::abs(cpu.startNll) &&
                      std::abs(cuda.minNll) < 1e-3 * std::abs(cuda.startNll);
  return Report(
      "fits with offsetting converge and agree (" + label + ")",
      offsetUsed && converged && worst <= kMaxFitSigmaDifference,
      TString::Format("CPU / CUDA: status %d / %d, EDM %.2g / %.2g (limit "
                      "%.0e), minimum %.4f / %.4f of an absolute NLL of "
                      "%.6g; max parameter difference %.3g sigma (limit %.2f)",
                      cpu.status, cuda.status, cpu.edm, cuda.edm,
                      kMaxConvergedEdm, cpu.minNll, cuda.minNll, cpu.startNll,
                      worst, kMaxFitSigmaDifference));
}

} // namespace

int main() {
#if !defined(AU_ROOFIT_BACKEND_CUDA) || !AU_ROOFIT_BACKEND_CUDA
  std::cout << "cuda-regression must be built with -DAU_ROOFIT_BACKEND_CUDA=1"
            << std::endl;
  return 2;
#else
  // As RooFitUtils does. ROOT master reports RooArgSet::checkForDup for the
  // range-projection nodes of a fixCoefRange'd simultaneous model; the direct
  // likelihood check below shows they do not affect the NLL.
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  for (Int_t i = 0; i < RooMsgService::instance().numStreams(); i++)
    RooMsgService::instance().getStream(i).removeTopic(RooFit::InputArguments);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  std::cout << "CUDA backend regression check, ROOT " << gROOT->GetVersion()
            << std::endl;

  Bool_t passed = kTRUE;
  ToyModel single = BuildModel(1, kEventsPerChannel);
  passed &= CheckConfiguration(
      single, TString::Format("single, %lld events", kEventsPerChannel), kFALSE,
      kReferenceSingleStartNll, kTRUE);

  ToyModel simultaneous = BuildModel(kSimultaneousChannels, kEventsPerChannel);
  SetYieldRanges(simultaneous, kTRUE);
  passed &= CheckConfiguration(
      simultaneous,
      TString::Format("simultaneous with SplitRange, %d x %lld events",
                      kSimultaneousChannels, kEventsPerChannel),
      kTRUE, kReferenceSplitStartNll, kTRUE);
  SetYieldRanges(simultaneous, kFALSE);
  passed &= CheckConfiguration(
      simultaneous,
      TString::Format("simultaneous over one shared range, %d x %lld events",
                      kSimultaneousChannels, kEventsPerChannel),
      kFALSE, kReferenceSharedStartNll, kFALSE);

  std::cout << TString::Format("fits with offsetting, %lld events per channel",
                               kFitEventsPerChannel)
            << std::endl;
  ToyModel fitSingle = BuildFitModel(1, kFitEventsPerChannel);
  passed &= CheckOffsetFits(fitSingle, "single", kFALSE);
  ToyModel fitSimultaneous =
      BuildFitModel(kSimultaneousChannels, kFitEventsPerChannel);
  // No shared-range fit: each channel's data only covers its own window, so a
  // model normalised over the union of windows does not describe it and the
  // fit has flat directions. The offset machinery is the same in every range
  // configuration, so these two fits already exercise it.
  SetYieldRanges(fitSimultaneous, kTRUE);
  passed &= CheckOffsetFits(fitSimultaneous, "SplitRange", kTRUE);

  std::cout << (passed ? "All checks passed." : "Some checks FAILED.")
            << std::endl;
  return passed ? 0 : 1;
#endif
}
