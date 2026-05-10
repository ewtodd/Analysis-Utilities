#include "RooFitPhotopeakPdfs.hpp"

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <cmath>

ClassImp(RooStepShelf);
ClassImp(RooLowExpTail);
ClassImp(RooLowLinTail);
ClassImp(RooHighExpTail);

namespace {

Double_t SoftPlusNeg(Double_t z) {
  if (z > 0)
    return std::log1p(std::exp(-z));
  return -z + std::log1p(std::exp(z));
}

Double_t SigmoidNeg(Double_t z) {
  if (z < 0)
    return 1.0 / (1.0 + std::exp(z));
  return std::exp(-z) / (1.0 + std::exp(-z));
}

Double_t StepAntideriv(Double_t z) {
  return -SoftPlusNeg(z) + SigmoidNeg(z);
}

Double_t ExpTailDensity(Double_t y, Double_t sigma, Double_t tau) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  return std::exp(y / tau) * std::erfc(y / sqrt2_sigma);
}

Double_t ExpTailAntideriv(Double_t y, Double_t sigma, Double_t tau) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t offset = sigma * sigma / tau;
  Double_t gauss_corr = std::exp(sigma * sigma / (2.0 * tau * tau));
  return tau * (std::exp(y / tau) * std::erfc(y / sqrt2_sigma) +
                gauss_corr * std::erf((y - offset) / sqrt2_sigma));
}

Double_t LowLinDensity(Double_t y, Double_t sigma, Double_t slope) {
  Double_t lin = 1.0 + slope * y;
  if (lin < 0.0)
    lin = 0.0;
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  return lin * std::erfc(y / sqrt2_sigma);
}

Double_t LowLinAntideriv(Double_t y, Double_t sigma, Double_t slope) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t beta = 1.0 / sqrt2_sigma;
  Double_t by = beta * y;
  Double_t erfc_by = std::erfc(by);
  Double_t erf_by = std::erf(by);
  Double_t gauss = std::exp(-(by * by));
  Double_t inv_beta_sqrt_pi = 1.0 / (beta * std::sqrt(M_PI));
  Double_t flat = y * erfc_by - gauss * inv_beta_sqrt_pi;
  Double_t lin = 0.5 * y * y * erfc_by -
                 y * gauss * inv_beta_sqrt_pi / 2.0 +
                 erf_by / (4.0 * beta * beta);
  return flat + slope * lin;
}

Double_t LowLinIntegral(Double_t y_lo, Double_t y_hi, Double_t sigma,
                        Double_t slope) {
  Double_t eff_lo = y_lo;
  Double_t eff_hi = y_hi;
  if (slope > 1e-10) {
    Double_t threshold = -1.0 / slope;
    if (threshold > eff_lo)
      eff_lo = threshold;
    if (eff_lo >= eff_hi)
      return 0.0;
  } else if (slope < -1e-10) {
    Double_t threshold = -1.0 / slope;
    if (threshold < eff_hi)
      eff_hi = threshold;
    if (eff_lo >= eff_hi)
      return 0.0;
  }
  return LowLinAntideriv(eff_hi, sigma, slope) -
         LowLinAntideriv(eff_lo, sigma, slope);
}

} // namespace

RooStepShelf::RooStepShelf(const char *name, const char *title, RooAbsReal &x,
                            RooAbsReal &mu, RooAbsReal &sigma)
    : RooAbsPdf(name, title), x_("x", "x", this, x),
      mu_("mu", "mu", this, mu), sigma_("sigma", "sigma", this, sigma) {}

RooStepShelf::RooStepShelf(const RooStepShelf &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_) {}

Double_t RooStepShelf::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t z = ((Double_t)x_ - (Double_t)mu_) / sigma;
  Double_t s = SigmoidNeg(z);
  return s * s;
}

Int_t RooStepShelf::getAnalyticalIntegral(RooArgSet &allVars,
                                          RooArgSet &analVars,
                                          const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooStepShelf::analyticalIntegral(Int_t code,
                                          const char *rangeName) const {
  if (code != 1)
    return 0.0;
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t z_lo = (x_lo - mu) / sigma;
  Double_t z_hi = (x_hi - mu) / sigma;
  return sigma * (StepAntideriv(z_hi) - StepAntideriv(z_lo));
}

RooLowExpTail::RooLowExpTail(const char *name, const char *title,
                              RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigma,
                              RooAbsReal &tau)
    : RooAbsPdf(name, title), x_("x", "x", this, x),
      mu_("mu", "mu", this, mu), sigma_("sigma", "sigma", this, sigma),
      tau_("tau", "tau", this, tau) {}

RooLowExpTail::RooLowExpTail(const RooLowExpTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      tau_("tau", this, other.tau_) {}

Double_t RooLowExpTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  Double_t tau = (Double_t)tau_;
  if (sigma <= 0 || tau <= 0)
    return 0.0;
  Double_t y = (Double_t)x_ - (Double_t)mu_;
  return ExpTailDensity(y, sigma, tau);
}

Int_t RooLowExpTail::getAnalyticalIntegral(RooArgSet &allVars,
                                            RooArgSet &analVars,
                                            const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooLowExpTail::analyticalIntegral(Int_t code,
                                            const char *rangeName) const {
  if (code != 1)
    return 0.0;
  Double_t sigma = (Double_t)sigma_;
  Double_t tau = (Double_t)tau_;
  if (sigma <= 0 || tau <= 0)
    return 0.0;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t y_lo = x_lo - mu;
  Double_t y_hi = x_hi - mu;
  return ExpTailAntideriv(y_hi, sigma, tau) -
         ExpTailAntideriv(y_lo, sigma, tau);
}

RooLowLinTail::RooLowLinTail(const char *name, const char *title,
                              RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigma,
                              RooAbsReal &slope)
    : RooAbsPdf(name, title), x_("x", "x", this, x),
      mu_("mu", "mu", this, mu), sigma_("sigma", "sigma", this, sigma),
      slope_("slope", "slope", this, slope) {}

RooLowLinTail::RooLowLinTail(const RooLowLinTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      slope_("slope", this, other.slope_) {}

Double_t RooLowLinTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t y = (Double_t)x_ - (Double_t)mu_;
  Double_t slope = (Double_t)slope_;
  return LowLinDensity(y, sigma, slope);
}

Int_t RooLowLinTail::getAnalyticalIntegral(RooArgSet &allVars,
                                            RooArgSet &analVars,
                                            const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooLowLinTail::analyticalIntegral(Int_t code,
                                            const char *rangeName) const {
  if (code != 1)
    return 0.0;
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t mu = (Double_t)mu_;
  Double_t slope = (Double_t)slope_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t y_lo = x_lo - mu;
  Double_t y_hi = x_hi - mu;
  return LowLinIntegral(y_lo, y_hi, sigma, slope);
}

RooHighExpTail::RooHighExpTail(const char *name, const char *title,
                                RooAbsReal &x, RooAbsReal &mu,
                                RooAbsReal &sigma, RooAbsReal &tau)
    : RooAbsPdf(name, title), x_("x", "x", this, x),
      mu_("mu", "mu", this, mu), sigma_("sigma", "sigma", this, sigma),
      tau_("tau", "tau", this, tau) {}

RooHighExpTail::RooHighExpTail(const RooHighExpTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      tau_("tau", this, other.tau_) {}

Double_t RooHighExpTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  Double_t tau = (Double_t)tau_;
  if (sigma <= 0 || tau <= 0)
    return 0.0;
  Double_t z = (Double_t)mu_ - (Double_t)x_;
  return ExpTailDensity(z, sigma, tau);
}

Int_t RooHighExpTail::getAnalyticalIntegral(RooArgSet &allVars,
                                             RooArgSet &analVars,
                                             const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooHighExpTail::analyticalIntegral(Int_t code,
                                             const char *rangeName) const {
  if (code != 1)
    return 0.0;
  Double_t sigma = (Double_t)sigma_;
  Double_t tau = (Double_t)tau_;
  if (sigma <= 0 || tau <= 0)
    return 0.0;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t z_lo = mu - x_lo;
  Double_t z_hi = mu - x_hi;
  return ExpTailAntideriv(z_lo, sigma, tau) -
         ExpTailAntideriv(z_hi, sigma, tau);
}
