import pickle
from functools import partial
from pathlib import Path

import numpy as np
from scipy.special import erf, erfc
from iminuit import Minuit, cost

import ROOT

_SQRT2 = np.sqrt(2.0)
_SQRT_2_OVER_PI = np.sqrt(2.0 / np.pi)
_SQRT_PI_OVER_2 = np.sqrt(np.pi / 2.0)
_INV_SQRT_2PI = 1.0 / np.sqrt(2.0 * np.pi)
_MAX_SUBSIDIARY_FRAC = 0.25

# ---------------------------------------------------------------------------------
# Analytic integrals for normalization computed using integral calculator the goat
# ---------------------------------------------------------------------------------


def _integral_gaussian(a, b, mu, sigma, amplitude):
    """Integral of amplitude * exp(-0.5*((x-mu)/sigma)^2) over [a, b]."""
    if amplitude == 0.0 or sigma <= 0.0:
        return 0.0
    return amplitude * sigma * _SQRT_PI_OVER_2 * (erf(
        (b - mu) / (_SQRT2 * sigma)) - erf((a - mu) / (_SQRT2 * sigma)))


def _integral_linear_background(a, b, bkg_constant, lin_bkg_slope):
    """Integral of (lin_bkg_slope * x + bkg_constant) over [a, b]."""
    return bkg_constant * (b - a) + lin_bkg_slope * (b * b - a * a) / 2.0


def _step_antideriv(u):
    """Antiderivative of 1/(1+exp(u))^2.
 
    F(u) = u - ln(1 + exp(u)) + 1/(1 + exp(u))
    """
    log1pexp = np.logaddexp(0.0, u)
    u_safe = np.clip(u, -500.0, 500.0)
    return u - log1pexp + 1.0 / (1.0 + np.exp(u_safe))


def _integral_step(a, b, mu, sigma, step_amplitude):
    """Integral of step_amplitude / (1+exp((x-mu)/sigma))^2 over [a, b]."""
    if step_amplitude == 0.0 or sigma <= 0.0:
        return 0.0
    za = (a - mu) / sigma
    zb = (b - mu) / sigma
    return step_amplitude * sigma * (_step_antideriv(zb) - _step_antideriv(za))


def _exp_erfc_antideriv(y, sigma, tau):
    """Antiderivative of exp(y/tau) * erfc(y/(sqrt(2)*sigma)).
 
    F(y) = tau * exp(sigma^2/(2*tau^2)) * erf(y/(sqrt(2)*sigma) - sigma/(sqrt(2)*tau))
         + tau * exp(y/tau) * erfc(y/(sqrt(2)*sigma))
    """
    u = y / (_SQRT2 * sigma)
    alpha = sigma / (_SQRT2 * tau)
    exp_arg = np.clip(alpha * alpha, -700.0, 700.0)
    term1 = tau * np.exp(exp_arg) * erf(u - alpha)
    term2 = tau * np.exp(np.clip(y / tau, -700.0, 700.0)) * erfc(u)
    return term1 + term2


def _erfc_const_antideriv(y, sigma):
    """Antiderivative of erfc(y/(sqrt(2)*sigma)).
 
    F(y) = y * erfc(y/(sqrt(2)*sigma)) - sigma*sqrt(2/pi)*exp(-y^2/(2*sigma^2))
    """
    u = y / (_SQRT2 * sigma)
    return y * erfc(u) - sigma * _SQRT_2_OVER_PI * np.exp(
        np.clip(-u * u, -700.0, 0.0))


def _erfc_slope_antideriv(y, sigma):
    """Antiderivative of y * erfc(y/(sqrt(2)*sigma)).
 
    F(y) = sigma^2/2 * erf(y/(sqrt(2)*sigma))
         + y^2/2 * erfc(y/(sqrt(2)*sigma))
         - y*sigma/sqrt(2*pi) * exp(-y^2/(2*sigma^2))
    """
    u = y / (_SQRT2 * sigma)
    s2 = sigma * sigma
    exp_term = np.exp(np.clip(-u * u, -700.0, 0.0))
    return (s2 / 2.0 * erf(u) + y * y / 2.0 * erfc(u) -
            y * sigma * _INV_SQRT_2PI * exp_term)


def _integral_low_tail(a, b, mu, sigma, exp_amp, exp_decay, lin_amp,
                       lin_slope):
    """Integral of low tail over [a, b] (y = x - mu)."""
    if sigma <= 0.0 or (exp_amp == 0.0 and lin_amp == 0.0):
        return 0.0
    ya = a - mu
    yb = b - mu
    result = 0.0
    if exp_amp != 0.0:
        result += exp_amp * (_exp_erfc_antideriv(yb, sigma, exp_decay) -
                             _exp_erfc_antideriv(ya, sigma, exp_decay))
    if lin_amp != 0.0:
        result += lin_amp * (_erfc_const_antideriv(yb, sigma) -
                             _erfc_const_antideriv(ya, sigma))
        result += lin_amp * lin_slope * (_erfc_slope_antideriv(yb, sigma) -
                                         _erfc_slope_antideriv(ya, sigma))
    return result


def _integral_high_tail(a, b, mu, sigma, exp_amp, exp_decay):
    """Integral of high tail over [a, b] (y = mu - x)."""
    if sigma <= 0.0 or exp_amp == 0.0:
        return 0.0
    ya = mu - b
    yb = mu - a
    return exp_amp * (_exp_erfc_antideriv(yb, sigma, exp_decay) -
                      _exp_erfc_antideriv(ya, sigma, exp_decay))


def _analytic_norm_peak(a, b, mu, sigma, gaus_amp, step_amp, low_exp_amp,
                        low_exp_decay, low_lin_amp, low_lin_slope,
                        high_exp_amp, high_exp_decay):
    """Analytic integral of one peak's components (no background) over [a, b]."""
    return (_integral_gaussian(a, b, mu, sigma, gaus_amp) +
            _integral_step(a, b, mu, sigma, step_amp) +
            _integral_low_tail(a, b, mu, sigma, low_exp_amp, low_exp_decay,
                               low_lin_amp, low_lin_slope) +
            _integral_high_tail(a, b, mu, sigma, high_exp_amp, high_exp_decay))


def _analytic_norm_single(a, b, mu, sigma, gaus_amp, step_amp, low_exp_amp,
                          low_exp_decay, low_lin_amp, low_lin_slope,
                          high_exp_amp, high_exp_decay, bkg_constant,
                          lin_bkg_slope):
    """Analytic integral of the full single-peak model over [a, b]."""
    return (_analytic_norm_peak(a, b, mu, sigma, gaus_amp, step_amp,
                                low_exp_amp, low_exp_decay, low_lin_amp,
                                low_lin_slope, high_exp_amp, high_exp_decay) +
            _integral_linear_background(a, b, bkg_constant, lin_bkg_slope))


def _analytic_norm_double(a, b, mu1, sigma1, gaus_amp1, step_amp1,
                          low_exp_amp1, low_exp_decay1, low_lin_amp1,
                          low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
                          sigma2, gaus_amp2, step_amp2, low_exp_amp2,
                          low_exp_decay2, low_lin_amp2, low_lin_slope2,
                          high_exp_amp2, high_exp_decay2, bkg_constant,
                          lin_bkg_slope):
    """Analytic integral of the full double-peak model over [a, b]."""
    return (
        _analytic_norm_peak(a, b, mu1, sigma1, gaus_amp1, step_amp1,
                            low_exp_amp1, low_exp_decay1, low_lin_amp1,
                            low_lin_slope1, high_exp_amp1, high_exp_decay1) +
        _analytic_norm_peak(a, b, mu2, sigma2, gaus_amp2, step_amp2,
                            low_exp_amp2, low_exp_decay2, low_lin_amp2,
                            low_lin_slope2, high_exp_amp2, high_exp_decay2) +
        _integral_linear_background(a, b, bkg_constant, lin_bkg_slope))


def _analytic_norm_triple(a, b, mu1, sigma1, gaus_amp1, step_amp1,
                          low_exp_amp1, low_exp_decay1, low_lin_amp1,
                          low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
                          sigma2, gaus_amp2, step_amp2, low_exp_amp2,
                          low_exp_decay2, low_lin_amp2, low_lin_slope2,
                          high_exp_amp2, high_exp_decay2, mu3, sigma3,
                          gaus_amp3, step_amp3, low_exp_amp3, low_exp_decay3,
                          low_lin_amp3, low_lin_slope3, high_exp_amp3,
                          high_exp_decay3, bkg_constant, lin_bkg_slope):
    """Analytic integral of the full triple-peak model over [a, b]."""
    return (
        _analytic_norm_peak(a, b, mu1, sigma1, gaus_amp1, step_amp1,
                            low_exp_amp1, low_exp_decay1, low_lin_amp1,
                            low_lin_slope1, high_exp_amp1, high_exp_decay1) +
        _analytic_norm_peak(a, b, mu2, sigma2, gaus_amp2, step_amp2,
                            low_exp_amp2, low_exp_decay2, low_lin_amp2,
                            low_lin_slope2, high_exp_amp2, high_exp_decay2) +
        _analytic_norm_peak(a, b, mu3, sigma3, gaus_amp3, step_amp3,
                            low_exp_amp3, low_exp_decay3, low_lin_amp3,
                            low_lin_slope3, high_exp_amp3, high_exp_decay3) +
        _integral_linear_background(a, b, bkg_constant, lin_bkg_slope))


_SINGLE_PEAK_PARAM_ORDER = [
    "mu", "sigma", "gaus_amp", "step_frac", "low_exp_frac", "low_exp_decay",
    "low_lin_frac", "low_lin_slope", "high_exp_frac", "high_exp_decay",
    "bkg_constant", "lin_bkg_slope"
]

_DOUBLE_PEAK_PARAM_ORDER = [
    "mu1", "sigma1", "gaus_amp1", "step_frac1", "low_exp_frac1",
    "low_exp_decay1", "low_lin_frac1", "low_lin_slope1", "high_exp_frac1",
    "high_exp_decay1", "mu2", "sigma2", "gaus_amp2", "step_frac2",
    "low_exp_frac2", "low_exp_decay2", "low_lin_frac2", "low_lin_slope2",
    "high_exp_frac2", "high_exp_decay2", "bkg_constant", "lin_bkg_slope"
]


def _gaussian(x, mu, sigma, amplitude):
    """Gaussian peak. Mirrors FittingFunctions::Gaussian."""
    z = (x - mu) / sigma
    return amplitude * np.exp(-0.5 * z * z)


def _linear_background(x, bkg_constant, lin_bkg_slope):
    """Linear background. Mirrors FittingFunctions::LinearBackground."""
    return lin_bkg_slope * x + bkg_constant


def _step(x, mu, sigma, step_amplitude):
    """Logistic step function. Mirrors FittingFunctions::Step."""
    if sigma <= 0.0:
        return np.zeros_like(x, dtype=np.float64)
    z = np.clip((x - mu) / sigma, -500, 350)
    denom = (1.0 + np.exp(z))**2
    return step_amplitude / denom


def _low_tail(x, mu, sigma, exp_amp, exp_decay, lin_amp, lin_slope):
    """Low-energy tail (exponential + linear) * erfc.

    Mirrors FittingFunctions::LowTail.
    """
    if sigma <= 0.0 or (exp_amp == 0.0 and lin_amp == 0.0):
        return np.zeros_like(x, dtype=np.float64)

    y = x - mu
    if exp_amp != 0.0:
        exp_term = exp_amp * np.exp(np.clip(y / exp_decay, -500, 500))
    else:
        exp_term = 0.0
    lin_term = lin_amp * (1.0 + lin_slope * y) if lin_amp != 0.0 else 0.0
    erfc_term = erfc(y / (_SQRT2 * sigma))

    return (exp_term + lin_term) * erfc_term


def _high_tail(x, mu, sigma, exp_amp, exp_decay):
    """High-energy exponential tail * erfc.

    Mirrors FittingFunctions::HighTail.
    """
    if sigma <= 0.0 or exp_amp == 0.0:
        return np.zeros_like(x, dtype=np.float64)

    y = mu - x
    return exp_amp * np.exp(np.clip(y / exp_decay, -500, 500)) * erfc(
        y / (_SQRT2 * sigma))


def _peak_function(x, mu, sigma, gaus_amp, step_amp, low_exp_amp,
                   low_exp_decay, low_lin_amp, low_lin_slope, high_exp_amp,
                   high_exp_decay, bkg_constant, lin_bkg_slope):
    """Full single-peak model (12 parameters).

    Mirrors FittingFunctions::PeakFunction.
    """
    return (_gaussian(x, mu, sigma, gaus_amp) +
            _linear_background(x, bkg_constant, lin_bkg_slope) +
            _step(x, mu, sigma, step_amp) +
            _low_tail(x, mu, sigma, low_exp_amp, low_exp_decay, low_lin_amp,
                      low_lin_slope) +
            _high_tail(x, mu, sigma, high_exp_amp, high_exp_decay))


def _double_peak_function(x, mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1,
                          low_exp_decay1, low_lin_amp1, low_lin_slope1,
                          high_exp_amp1, high_exp_decay1, mu2, sigma2,
                          gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
                          low_lin_amp2, low_lin_slope2, high_exp_amp2,
                          high_exp_decay2, bkg_constant, lin_bkg_slope):
    """Full double-peak model (22 parameters, shared background).

    Mirrors FittingFunctions::DoublePeakFunction.
    Peak 1: params 0-9, Peak 2: params 10-19, Background: params 20-21.
    """
    return (_gaussian(x, mu1, sigma1, gaus_amp1) +
            _step(x, mu1, sigma1, step_amp1) +
            _low_tail(x, mu1, sigma1, low_exp_amp1, low_exp_decay1,
                      low_lin_amp1, low_lin_slope1) +
            _high_tail(x, mu1, sigma1, high_exp_amp1, high_exp_decay1) +
            _gaussian(x, mu2, sigma2, gaus_amp2) +
            _step(x, mu2, sigma2, step_amp2) +
            _low_tail(x, mu2, sigma2, low_exp_amp2, low_exp_decay2,
                      low_lin_amp2, low_lin_slope2) +
            _high_tail(x, mu2, sigma2, high_exp_amp2, high_exp_decay2) +
            _linear_background(x, bkg_constant, lin_bkg_slope))


def _triple_peak_function(x, mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1,
                          low_exp_decay1, low_lin_amp1, low_lin_slope1,
                          high_exp_amp1, high_exp_decay1, mu2, sigma2,
                          gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
                          low_lin_amp2, low_lin_slope2, high_exp_amp2,
                          high_exp_decay2, mu3, sigma3, gaus_amp3, step_amp3,
                          low_exp_amp3, low_exp_decay3, low_lin_amp3,
                          low_lin_slope3, high_exp_amp3, high_exp_decay3,
                          bkg_constant, lin_bkg_slope):
    """Full triple-peak model (32 parameters, shared background).

    Mirrors FittingFunctions::TriplePeakFunction.
    Peak 1: 0-9, Peak 2: 10-19, Peak 3: 20-29, Background: 30-31.
    """
    return (_gaussian(x, mu1, sigma1, gaus_amp1) +
            _step(x, mu1, sigma1, step_amp1) +
            _low_tail(x, mu1, sigma1, low_exp_amp1, low_exp_decay1,
                      low_lin_amp1, low_lin_slope1) +
            _high_tail(x, mu1, sigma1, high_exp_amp1, high_exp_decay1) +
            _gaussian(x, mu2, sigma2, gaus_amp2) +
            _step(x, mu2, sigma2, step_amp2) +
            _low_tail(x, mu2, sigma2, low_exp_amp2, low_exp_decay2,
                      low_lin_amp2, low_lin_slope2) +
            _high_tail(x, mu2, sigma2, high_exp_amp2, high_exp_decay2) +
            _gaussian(x, mu3, sigma3, gaus_amp3) +
            _step(x, mu3, sigma3, step_amp3) +
            _low_tail(x, mu3, sigma3, low_exp_amp3, low_exp_decay3,
                      low_lin_amp3, low_lin_slope3) +
            _high_tail(x, mu3, sigma3, high_exp_amp3, high_exp_decay3) +
            _linear_background(x, bkg_constant, lin_bkg_slope))


def _estimate_from_data(data, fit_range_low, fit_range_high):
    """Histogram unbinned data and extract peak height + background estimate.

    Mirrors the C++ constructor logic and EstimateBackground().
    """
    in_range = data[(data >= fit_range_low) & (data <= fit_range_high)]
    n_bins = max(10, int(np.sqrt(len(in_range))))
    counts, edges = np.histogram(in_range,
                                 bins=n_bins,
                                 range=(fit_range_low, fit_range_high))

    peak_height = float(counts.max())
    n_sideband = max(1, n_bins // 10)
    left_avg = counts[:n_sideband].mean()
    right_avg = counts[-n_sideband:].mean()
    bkg_estimate = float((left_avg + right_avg) / 2.0)

    return peak_height, bkg_estimate


def _estimate_single_peak(data,
                          fit_range_low,
                          fit_range_high,
                          peak_height,
                          bkg_estimate,
                          use_step,
                          use_low_exp_tail,
                          use_low_lin_tail,
                          use_high_exp_tail,
                          use_flat_background,
                          suffix=""):
    """Build values/limits/fixed dicts for one peak plus background.

    suffix is "" for single peak, "1"/"2"/"3" for multi-peak.
    """
    range_width = fit_range_high - fit_range_low

    values = {}
    limits = {}
    fixed = {}

    # Gaussian
    values[f"mu{suffix}"] = (fit_range_low + fit_range_high) / 2.0
    limits[f"mu{suffix}"] = (fit_range_low, fit_range_high)
    values[f"sigma{suffix}"] = range_width * 0.01
    limits[f"sigma{suffix}"] = (range_width * 0.001, range_width * 0.5)
    values[f"gaus_amp{suffix}"] = peak_height * 0.999
    limits[f"gaus_amp{suffix}"] = (0, peak_height * 2.0)

    # Step fraction
    if use_step:
        values[f"step_frac{suffix}"] = 0.1
        limits[f"step_frac{suffix}"] = (0, _MAX_SUBSIDIARY_FRAC)
    else:
        values[f"step_frac{suffix}"] = 0.0
        fixed[f"step_frac{suffix}"] = True

    # Low exponential tail fraction
    if use_low_exp_tail:
        values[f"low_exp_frac{suffix}"] = 0.15
        limits[f"low_exp_frac{suffix}"] = (0, _MAX_SUBSIDIARY_FRAC)
        values[f"low_exp_decay{suffix}"] = 1.0
        limits[f"low_exp_decay{suffix}"] = (0.1, 50)
    else:
        values[f"low_exp_frac{suffix}"] = 0.0
        values[f"low_exp_decay{suffix}"] = 1.0
        fixed[f"low_exp_frac{suffix}"] = True
        fixed[f"low_exp_decay{suffix}"] = True

    # Low linear tail fraction
    if use_low_lin_tail:
        values[f"low_lin_frac{suffix}"] = 0.15
        limits[f"low_lin_frac{suffix}"] = (0, _MAX_SUBSIDIARY_FRAC)
        values[f"low_lin_slope{suffix}"] = 0.0
        limits[f"low_lin_slope{suffix}"] = (-0.5 * bkg_estimate / range_width,
                                            0.5 * bkg_estimate / range_width)
    else:
        values[f"low_lin_frac{suffix}"] = 0.0
        values[f"low_lin_slope{suffix}"] = 0.0
        fixed[f"low_lin_frac{suffix}"] = True
        fixed[f"low_lin_slope{suffix}"] = True

    # High exponential tail fraction
    if use_high_exp_tail:
        values[f"high_exp_frac{suffix}"] = 0.15
        limits[f"high_exp_frac{suffix}"] = (0, _MAX_SUBSIDIARY_FRAC)
        values[f"high_exp_decay{suffix}"] = 1.0
        limits[f"high_exp_decay{suffix}"] = (0.1, 50)
    else:
        values[f"high_exp_frac{suffix}"] = 0.0
        values[f"high_exp_decay{suffix}"] = 1.0
        fixed[f"high_exp_frac{suffix}"] = True
        fixed[f"high_exp_decay{suffix}"] = True

    return values, limits, fixed


def estimate_peak_params(data,
                         fit_range_low,
                         fit_range_high,
                         use_step=False,
                         use_low_exp_tail=False,
                         use_low_lin_tail=False,
                         use_high_exp_tail=False,
                         use_flat_background=False,
                         **overrides):
    """Estimate initial parameters for a single-peak fit.

    Returns (values, limits, fixed) dicts for use with iminuit::

        values, limits, fixed = estimate_peak_params(data, 640, 680, mu=661.7)
        m = Minuit(c, **values)
        for k, v in limits.items(): m.limits[k] = v
        for k, v in fixed.items(): m.fixed[k] = v

    Any keyword argument overrides the corresponding auto-estimated value.
    """
    peak_height, bkg_estimate = _estimate_from_data(data, fit_range_low,
                                                    fit_range_high)
    range_width = fit_range_high - fit_range_low

    values, limits, fixed = _estimate_single_peak(
        data, fit_range_low, fit_range_high, peak_height, bkg_estimate,
        use_step, use_low_exp_tail, use_low_lin_tail, use_high_exp_tail,
        use_flat_background)

    # Background (shared, no suffix)
    values["bkg_constant"] = bkg_estimate
    limits["bkg_constant"] = (0, peak_height * 2.0)

    if use_flat_background:
        values["lin_bkg_slope"] = 0.0
        fixed["lin_bkg_slope"] = True
    else:
        values["lin_bkg_slope"] = 0.0
        limits["lin_bkg_slope"] = (-0.5 * bkg_estimate / range_width,
                                   0.5 * bkg_estimate / range_width)

    values.update(overrides)
    return values, limits, fixed


def estimate_double_peak_params(data,
                                fit_range_low,
                                fit_range_high,
                                use_step=False,
                                use_low_exp_tail=False,
                                use_low_lin_tail=False,
                                use_high_exp_tail=False,
                                use_flat_background=False,
                                **overrides):
    """Estimate initial parameters for a double-peak fit.

    You almost certainly want to override mu1 and mu2.
    Returns (values, limits, fixed) dicts for use with iminuit.
    """
    peak_height, bkg_estimate = _estimate_from_data(data, fit_range_low,
                                                    fit_range_high)
    range_width = fit_range_high - fit_range_low

    v1, l1, f1 = _estimate_single_peak(data,
                                       fit_range_low,
                                       fit_range_high,
                                       peak_height,
                                       bkg_estimate,
                                       use_step,
                                       use_low_exp_tail,
                                       use_low_lin_tail,
                                       use_high_exp_tail,
                                       use_flat_background,
                                       suffix="1")
    v2, l2, f2 = _estimate_single_peak(data,
                                       fit_range_low,
                                       fit_range_high,
                                       peak_height,
                                       bkg_estimate,
                                       use_step,
                                       use_low_exp_tail,
                                       use_low_lin_tail,
                                       use_high_exp_tail,
                                       use_flat_background,
                                       suffix="2")

    values = {**v1, **v2}
    limits = {**l1, **l2}
    fixed = {**f1, **f2}

    values["bkg_constant"] = bkg_estimate
    limits["bkg_constant"] = (0, peak_height * 2.0)
    if use_flat_background:
        values["lin_bkg_slope"] = 0.0
        fixed["lin_bkg_slope"] = True
    else:
        values["lin_bkg_slope"] = 0.0
        limits["lin_bkg_slope"] = (-0.5 * bkg_estimate / range_width,
                                   0.5 * bkg_estimate / range_width)

    values.update(overrides)
    return values, limits, fixed


def estimate_triple_peak_params(data,
                                fit_range_low,
                                fit_range_high,
                                use_step=False,
                                use_low_exp_tail=False,
                                use_low_lin_tail=False,
                                use_high_exp_tail=False,
                                use_flat_background=False,
                                **overrides):
    """Estimate initial parameters for a triple-peak fit.

    You almost certainly want to override mu1, mu2, and mu3.
    Returns (values, limits, fixed) dicts for use with iminuit.
    """
    peak_height, bkg_estimate = _estimate_from_data(data, fit_range_low,
                                                    fit_range_high)
    range_width = fit_range_high - fit_range_low

    v1, l1, f1 = _estimate_single_peak(data,
                                       fit_range_low,
                                       fit_range_high,
                                       peak_height,
                                       bkg_estimate,
                                       use_step,
                                       use_low_exp_tail,
                                       use_low_lin_tail,
                                       use_high_exp_tail,
                                       use_flat_background,
                                       suffix="1")
    v2, l2, f2 = _estimate_single_peak(data,
                                       fit_range_low,
                                       fit_range_high,
                                       peak_height,
                                       bkg_estimate,
                                       use_step,
                                       use_low_exp_tail,
                                       use_low_lin_tail,
                                       use_high_exp_tail,
                                       use_flat_background,
                                       suffix="2")
    v3, l3, f3 = _estimate_single_peak(data,
                                       fit_range_low,
                                       fit_range_high,
                                       peak_height,
                                       bkg_estimate,
                                       use_step,
                                       use_low_exp_tail,
                                       use_low_lin_tail,
                                       use_high_exp_tail,
                                       use_flat_background,
                                       suffix="3")

    values = {**v1, **v2, **v3}
    limits = {**l1, **l2, **l3}
    fixed = {**f1, **f2, **f3}

    values["bkg_constant"] = bkg_estimate
    limits["bkg_constant"] = (0, peak_height * 2.0)
    if use_flat_background:
        values["lin_bkg_slope"] = 0.0
        fixed["lin_bkg_slope"] = True
    else:
        values["lin_bkg_slope"] = 0.0
        limits["lin_bkg_slope"] = (-0.5 * bkg_estimate / range_width,
                                   0.5 * bkg_estimate / range_width)

    values.update(overrides)
    return values, limits, fixed


def _single_peak_pdf(fit_range_low, fit_range_high, x, mu, sigma, gaus_amp,
                     step_frac, low_exp_frac, low_exp_decay, low_lin_frac,
                     low_lin_slope, high_exp_frac, high_exp_decay,
                     bkg_constant, lin_bkg_slope):
    step_amp = step_frac * gaus_amp
    low_exp_amp = low_exp_frac * gaus_amp
    low_lin_amp = low_lin_frac * gaus_amp
    high_exp_amp = high_exp_frac * gaus_amp
    args = (mu, sigma, gaus_amp, step_amp, low_exp_amp, low_exp_decay,
            low_lin_amp, low_lin_slope, high_exp_amp, high_exp_decay,
            bkg_constant, lin_bkg_slope)
    norm = _analytic_norm_single(fit_range_low, fit_range_high, *args)
    return _peak_function(x, *args) / norm


def _double_peak_pdf(fit_range_low, fit_range_high, x, mu1, sigma1, gaus_amp1,
                     step_frac1, low_exp_frac1, low_exp_decay1, low_lin_frac1,
                     low_lin_slope1, high_exp_frac1, high_exp_decay1, mu2,
                     sigma2, gaus_amp2, step_frac2, low_exp_frac2,
                     low_exp_decay2, low_lin_frac2, low_lin_slope2,
                     high_exp_frac2, high_exp_decay2, bkg_constant,
                     lin_bkg_slope):
    step_amp1 = step_frac1 * gaus_amp1
    low_exp_amp1 = low_exp_frac1 * gaus_amp1
    low_lin_amp1 = low_lin_frac1 * gaus_amp1
    high_exp_amp1 = high_exp_frac1 * gaus_amp1
    step_amp2 = step_frac2 * gaus_amp2
    low_exp_amp2 = low_exp_frac2 * gaus_amp2
    low_lin_amp2 = low_lin_frac2 * gaus_amp2
    high_exp_amp2 = high_exp_frac2 * gaus_amp2
    args = (mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1, low_exp_decay1,
            low_lin_amp1, low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
            sigma2, gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
            low_lin_amp2, low_lin_slope2, high_exp_amp2, high_exp_decay2,
            bkg_constant, lin_bkg_slope)
    norm = _analytic_norm_double(fit_range_low, fit_range_high, *args)
    return _double_peak_function(x, *args) / norm


def _triple_peak_pdf(fit_range_low, fit_range_high, x, mu1, sigma1, gaus_amp1,
                     step_frac1, low_exp_frac1, low_exp_decay1, low_lin_frac1,
                     low_lin_slope1, high_exp_frac1, high_exp_decay1, mu2,
                     sigma2, gaus_amp2, step_frac2, low_exp_frac2,
                     low_exp_decay2, low_lin_frac2, low_lin_slope2,
                     high_exp_frac2, high_exp_decay2, mu3, sigma3, gaus_amp3,
                     step_frac3, low_exp_frac3, low_exp_decay3, low_lin_frac3,
                     low_lin_slope3, high_exp_frac3, high_exp_decay3,
                     bkg_constant, lin_bkg_slope):
    step_amp1 = step_frac1 * gaus_amp1
    low_exp_amp1 = low_exp_frac1 * gaus_amp1
    low_lin_amp1 = low_lin_frac1 * gaus_amp1
    high_exp_amp1 = high_exp_frac1 * gaus_amp1
    step_amp2 = step_frac2 * gaus_amp2
    low_exp_amp2 = low_exp_frac2 * gaus_amp2
    low_lin_amp2 = low_lin_frac2 * gaus_amp2
    high_exp_amp2 = high_exp_frac2 * gaus_amp2
    step_amp3 = step_frac3 * gaus_amp3
    low_exp_amp3 = low_exp_frac3 * gaus_amp3
    low_lin_amp3 = low_lin_frac3 * gaus_amp3
    high_exp_amp3 = high_exp_frac3 * gaus_amp3
    args = (mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1, low_exp_decay1,
            low_lin_amp1, low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
            sigma2, gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
            low_lin_amp2, low_lin_slope2, high_exp_amp2, high_exp_decay2, mu3,
            sigma3, gaus_amp3, step_amp3, low_exp_amp3, low_exp_decay3,
            low_lin_amp3, low_lin_slope3, high_exp_amp3, high_exp_decay3,
            bkg_constant, lin_bkg_slope)
    norm = _analytic_norm_triple(fit_range_low, fit_range_high, *args)
    return _triple_peak_function(x, *args) / norm


def single_peak_pdf(fit_range_low, fit_range_high):
    """Return a normalized single-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_single_peak_pdf, fit_range_low, fit_range_high)


def double_peak_pdf(fit_range_low, fit_range_high):
    """Return a normalized double-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_double_peak_pdf, fit_range_low, fit_range_high)


def triple_peak_pdf(fit_range_low, fit_range_high):
    """Return a normalized triple-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_triple_peak_pdf, fit_range_low, fit_range_high)


def _frac_to_amp_params(minuit_result, param_order):
    """Convert fraction-parameterized Minuit values to absolute amplitudes.

    Returns a tuple of parameter values in the same order as param_order,
    but with *_frac entries replaced by frac * gaus_amp (the absolute amplitude
    that the C++ functions expect).
    """
    vals = []
    for name in param_order:
        v = minuit_result.values[name]
        if "_frac" in name:
            # Determine which peak's gaus_amp to multiply by
            suffix = name.split("_frac")[-1]  # "" or "1" or "2" or "3"
            gaus_key = f"gaus_amp{suffix}"
            v = v * minuit_result.values[gaus_key]
        vals.append(v)
    return tuple(vals)


def _frac_to_amp_errors(minuit_result, param_order):
    """Convert fraction-parameterized Minuit errors to absolute amplitude errors.

    Uses simple product rule: err(frac * gaus_amp) ≈ frac * gaus_amp_err
    + gaus_amp * frac_err (added in quadrature), but for plotting purposes
    the dominant term is gaus_amp * frac_err since gaus_amp is well-determined.
    """
    errs = []
    for name in param_order:
        e = minuit_result.errors[name]
        if "_frac" in name:
            suffix = name.split("_frac")[-1]
            gaus_key = f"gaus_amp{suffix}"
            gaus_amp = minuit_result.values[gaus_key]
            gaus_err = minuit_result.errors[gaus_key]
            frac = minuit_result.values[name]
            e = np.sqrt((gaus_amp * e)**2 + (frac * gaus_err)**2)
        errs.append(e)
    return tuple(errs)


# C++ parameter order (absolute amplitudes, not fractions)
_SINGLE_PEAK_CPP_ORDER = [
    "mu", "sigma", "gaus_amp", "step_amp", "low_exp_amp", "low_exp_decay",
    "low_lin_amp", "low_lin_slope", "high_exp_amp", "high_exp_decay",
    "bkg_constant", "lin_bkg_slope"
]

_DOUBLE_PEAK_CPP_ORDER = [
    "mu1", "sigma1", "gaus_amp1", "step_amp1", "low_exp_amp1",
    "low_exp_decay1", "low_lin_amp1", "low_lin_slope1", "high_exp_amp1",
    "high_exp_decay1", "mu2", "sigma2", "gaus_amp2", "step_amp2",
    "low_exp_amp2", "low_exp_decay2", "low_lin_amp2", "low_lin_slope2",
    "high_exp_amp2", "high_exp_decay2", "bkg_constant", "lin_bkg_slope"
]


def plot_single_peak_fit(hist, minuit_result, fit_range_low, fit_range_high,
                         input_name, peak_name):
    """Plot a single-peak fit result using the C++ PlotFitSinglePeak.

    Parameters
    ----------
    hist : ROOT.TH1
        The histogram that was fitted.
    minuit_result : iminuit.Minuit
        Completed Minuit object from fit_single_peak.
    fit_range_low, fit_range_high : float
        Fit range boundaries.
    input_name, peak_name : str
        Labels passed to PlotFitSinglePeak for the output filename.
    """
    # Compute scale factor: N_events * bin_width / integral(peak_function)
    # so the normalized PDF parameters match histogram counts.
    bin_low = hist.FindBin(fit_range_low)
    bin_high = hist.FindBin(fit_range_high)
    n_events = hist.Integral(bin_low, bin_high)
    bin_width = hist.GetBinWidth(1)

    # Convert fracs → absolute amps for norm computation and C++ TF1
    abs_params = _frac_to_amp_params(minuit_result, _SINGLE_PEAK_PARAM_ORDER)
    abs_errors = _frac_to_amp_errors(minuit_result, _SINGLE_PEAK_PARAM_ORDER)
    norm = _analytic_norm_single(fit_range_low, fit_range_high, *abs_params)
    scale = n_events * bin_width / norm

    # Indices of amplitude parameters (linear in the function value).
    # Shape params (mu, sigma, exp_decay, lin_slope) are unchanged.
    amp_indices = {2, 3, 4, 6, 8, 10, 11}

    def _active(frac_name):
        return (not minuit_result.fixed[frac_name]
                and abs(minuit_result.values[frac_name]) > 1e-6)

    fitter = ROOT.FittingUtils(hist, fit_range_low, fit_range_high, False,
                               _active("step_frac"), _active("low_exp_frac"),
                               _active("low_lin_frac"),
                               _active("high_exp_frac"))
    tf1 = fitter.GetFitFunction()
    for i in range(len(_SINGLE_PEAK_PARAM_ORDER)):
        val = abs_params[i]
        err = abs_errors[i]
        if i in amp_indices:
            val *= scale
            err *= scale
        tf1.SetParameter(i, val)
        tf1.SetParError(i, err)
        if minuit_result.fixed[_SINGLE_PEAK_PARAM_ORDER[i]]:
            tf1.FixParameter(i, val)

    # Chi2/ndf from binned data
    chi2 = 0.0
    n_bins_used = 0
    for i in range(bin_low, bin_high + 1):
        observed = hist.GetBinContent(i)
        x = hist.GetBinCenter(i)
        expected = tf1.Eval(x)
        if expected > 0:
            chi2 += (observed - expected)**2 / expected
            n_bins_used += 1
    n_free = sum(1 for name in _SINGLE_PEAK_PARAM_ORDER
                 if not minuit_result.fixed[name])
    ndf = n_bins_used - n_free
    if ndf > 0:
        chi2_ndf = chi2 / ndf
        chi2_label = f"#chi^{{2}}/ndf = {chi2_ndf:.3f}"
        print(f"Chi2/ndf = {chi2:.1f}/{ndf} = {chi2_ndf:.3f}")
    else:
        chi2_label = f"#chi^{{2}} = {chi2:.1f}"
        print(
            f"Chi2 = {chi2:.1f} (ndf <= 0, {n_bins_used} bins, {n_free} free params)"
        )

    fitter.PlotFitSinglePeak(input_name, peak_name, chi2_label)


def plot_double_peak_fit(hist, minuit_result, fit_range_low, fit_range_high,
                         input_name, peak_name):
    """Plot a double-peak fit result using the C++ PlotFitDoublePeak."""
    bin_low = hist.FindBin(fit_range_low)
    bin_high = hist.FindBin(fit_range_high)
    n_events = hist.Integral(bin_low, bin_high)
    bin_width = hist.GetBinWidth(1)

    # Convert fracs → absolute amps for norm computation and C++ TF1
    abs_params = _frac_to_amp_params(minuit_result, _DOUBLE_PEAK_PARAM_ORDER)
    abs_errors = _frac_to_amp_errors(minuit_result, _DOUBLE_PEAK_PARAM_ORDER)
    norm = _analytic_norm_double(fit_range_low, fit_range_high, *abs_params)
    scale = n_events * bin_width / norm

    amp_indices = {2, 3, 4, 6, 8, 12, 13, 14, 16, 18, 20, 21}

    def _active(frac_name):
        return (not minuit_result.fixed[frac_name]
                and abs(minuit_result.values[frac_name]) > 1e-6)

    fitter = ROOT.FittingUtils(
        hist, fit_range_low, fit_range_high, False,
        _active("step_frac1") or _active("step_frac2"),
        _active("low_exp_frac1") or _active("low_exp_frac2"),
        _active("low_lin_frac1") or _active("low_lin_frac2"),
        _active("high_exp_frac1") or _active("high_exp_frac2"))

    tf1 = ROOT.TF1("DoublePeak", ROOT.FittingFunctions.DoublePeakFunction,
                   fit_range_low, fit_range_high, 22)
    for i in range(len(_DOUBLE_PEAK_PARAM_ORDER)):
        val = abs_params[i]
        err = abs_errors[i]
        if i in amp_indices:
            val *= scale
            err *= scale
        tf1.SetParameter(i, val)
        tf1.SetParError(i, err)
        if minuit_result.fixed[_DOUBLE_PEAK_PARAM_ORDER[i]]:
            tf1.FixParameter(i, val)
    fitter.SetFitFunction(tf1)

    chi2 = 0.0
    n_bins_used = 0
    for i in range(bin_low, bin_high + 1):
        observed = hist.GetBinContent(i)
        x = hist.GetBinCenter(i)
        expected = tf1.Eval(x)
        if expected > 0:
            chi2 += (observed - expected)**2 / expected
            n_bins_used += 1
    n_free = sum(1 for name in _DOUBLE_PEAK_PARAM_ORDER
                 if not minuit_result.fixed[name])
    ndf = n_bins_used - n_free
    if ndf > 0:
        chi2_ndf = chi2 / ndf
        chi2_label = f"#chi^{{2}}/ndf = {chi2_ndf:.3f}"
        print(f"Chi2/ndf = {chi2:.1f}/{ndf} = {chi2_ndf:.3f}")
    else:
        chi2_label = f"#chi^{{2}} = {chi2:.1f}"
        print(
            f"Chi2 = {chi2:.1f} (ndf <= 0, {n_bins_used} bins, {n_free} free params)"
        )

    fitter.PlotFitDoublePeak(input_name, peak_name, chi2_label)


def _save_minuit_state(m, path):
    state = {
        "values": {
            k: m.values[k]
            for k in m.parameters
        },
        "errors": {
            k: m.errors[k]
            for k in m.parameters
        },
        "fixed": {
            k: m.fixed[k]
            for k in m.parameters
        },
        "limits": {
            k: m.limits[k]
            for k in m.parameters
        },
        "fval": m.fval,
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(state, f)
    print(f"Fit result cached to {path}")


def _restore_minuit_state(m, path):
    with open(path, "rb") as f:
        state = pickle.load(f)
    for k, v in state["values"].items():
        m.values[k] = v
    for k, v in state["errors"].items():
        m.errors[k] = v
    for k, v in state["limits"].items():
        m.limits[k] = v
    for k, v in state["fixed"].items():
        m.fixed[k] = v
    print(f"Loaded cached fit from {path} (NLL = {state['fval']:.2f})")


def fit_single_peak(df_column,
                    fit_range_low,
                    fit_range_high,
                    expected_mu,
                    cache_path=None):
    """Fit a single-peak model with automated component testing and pruning.

    Returns a Minuit object. If cache_path is given, results are cached.
    """
    data = np.asarray(df_column, dtype=np.float64)
    data = data[(data >= fit_range_low) & (data <= fit_range_high)]
    pdf = single_peak_pdf(fit_range_low, fit_range_high)
    c = cost.UnbinnedNLL(data, pdf)
    values, limits, fixed = estimate_peak_params(data,
                                                 fit_range_low,
                                                 fit_range_high,
                                                 mu=expected_mu)
    m = Minuit(c, **values)
    for k, v in limits.items():
        m.limits[k] = v
    for k, v in fixed.items():
        m.fixed[k] = v

    if cache_path and Path(cache_path).exists():
        _restore_minuit_state(m, cache_path)
        return m

    m.migrad()
    best_nll = m.fval
    gaus_amp = m.values["gaus_amp"]
    print(f"Initial fit: NLL = {best_nll:.2f}")

    print("Testing low-side group (step + low exp tail + low lin tail)...")
    m.fixed["step_frac"] = False
    m.limits["step_frac"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.values["step_frac"] = 0.1

    m.fixed["low_exp_frac"] = False
    m.fixed["low_exp_decay"] = False
    m.limits["low_exp_frac"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["low_exp_decay"] = (0.1, 50)
    m.values["low_exp_frac"] = 0.15
    m.values["low_exp_decay"] = 1.0

    m.fixed["low_lin_frac"] = False
    m.fixed["low_lin_slope"] = False
    m.limits["low_lin_frac"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["low_lin_slope"] = (-1, 1)
    m.values["low_lin_frac"] = 0.15
    m.values["low_lin_slope"] = 0.0

    m.migrad()
    print(
        f"  Low-side group NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
    )

    if m.fval < best_nll - 1:
        print("  Low-side group ACCEPTED")
        best_nll = m.fval
    else:
        print("  Low-side group REJECTED — re-fixing all")
        m.fixed["step_frac"] = True
        m.values["step_frac"] = 0.0
        m.fixed["low_exp_frac"] = True
        m.fixed["low_exp_decay"] = True
        m.values["low_exp_frac"] = 0.0
        m.values["low_exp_decay"] = 1.0
        m.fixed["low_lin_frac"] = True
        m.fixed["low_lin_slope"] = True
        m.values["low_lin_frac"] = 0.0
        m.values["low_lin_slope"] = 0.0

    print("Testing high exponential tail...")
    m.fixed["high_exp_frac"] = False
    m.fixed["high_exp_decay"] = False
    m.limits["high_exp_frac"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["high_exp_decay"] = (0.1, 50)
    m.values["high_exp_frac"] = 0.15
    m.values["high_exp_decay"] = 1.0
    m.migrad()
    print(f"  High tail NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})")
    if m.fval < best_nll - 1:
        best_nll = m.fval
        print("  High tail ACCEPTED")
    else:
        m.fixed["high_exp_frac"] = True
        m.fixed["high_exp_decay"] = True
        m.values["high_exp_frac"] = 0.0
        m.values["high_exp_decay"] = 1.0
        print("  High tail REJECTED")

    print("Final fit with selected components...")
    m.migrad()
    m.hesse()
    print(f"Final NLL = {m.fval:.2f}")
    print(m)
    if cache_path:
        _save_minuit_state(m, cache_path)
    return m


def _test_low_side_group(m, suffix, gaus_amp, best_nll):
    """Test and accept/reject low-side components for a single peak within a
    multi-peak fit. Returns updated best_nll."""
    s = suffix
    print(f"Testing low-side group for peak{s}...")

    m.fixed[f"step_frac{s}"] = False
    m.limits[f"step_frac{s}"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.values[f"step_frac{s}"] = 0.1

    m.fixed[f"low_exp_frac{s}"] = False
    m.fixed[f"low_exp_decay{s}"] = False
    m.limits[f"low_exp_frac{s}"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits[f"low_exp_decay{s}"] = (0.1, 50)
    m.values[f"low_exp_frac{s}"] = 0.15
    m.values[f"low_exp_decay{s}"] = 1.0

    m.fixed[f"low_lin_frac{s}"] = False
    m.fixed[f"low_lin_slope{s}"] = False
    m.limits[f"low_lin_frac{s}"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits[f"low_lin_slope{s}"] = (-1, 1)
    m.values[f"low_lin_frac{s}"] = 0.15
    m.values[f"low_lin_slope{s}"] = 0.0

    m.migrad()
    print(
        f"  Low-side group NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
    )

    if m.fval < best_nll - 1:
        print(f"  Low-side group peak{s} ACCEPTED")
        best_nll = m.fval
    else:
        print(f"  Low-side group peak{s} REJECTED — re-fixing all")
        m.fixed[f"step_frac{s}"] = True
        m.values[f"step_frac{s}"] = 0.0
        m.fixed[f"low_exp_frac{s}"] = True
        m.fixed[f"low_exp_decay{s}"] = True
        m.values[f"low_exp_frac{s}"] = 0.0
        m.values[f"low_exp_decay{s}"] = 1.0
        m.fixed[f"low_lin_frac{s}"] = True
        m.fixed[f"low_lin_slope{s}"] = True
        m.values[f"low_lin_frac{s}"] = 0.0
        m.values[f"low_lin_slope{s}"] = 0.0

    return best_nll


def _test_high_tail(m, suffix, gaus_amp, best_nll):
    """Test high exponential tail for a single peak. Returns updated best_nll."""
    s = suffix
    print(f"Testing high exponential tail for peak{s}...")
    m.fixed[f"high_exp_frac{s}"] = False
    m.fixed[f"high_exp_decay{s}"] = False
    m.limits[f"high_exp_frac{s}"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits[f"high_exp_decay{s}"] = (0.1, 50)
    m.values[f"high_exp_frac{s}"] = 0.15
    m.values[f"high_exp_decay{s}"] = 1.0
    m.migrad()
    print(
        f"  High tail{s} NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})")
    if m.fval < best_nll - 1:
        best_nll = m.fval
        print(f"  High tail{s} ACCEPTED")
    else:
        m.fixed[f"high_exp_frac{s}"] = True
        m.fixed[f"high_exp_decay{s}"] = True
        m.values[f"high_exp_frac{s}"] = 0.0
        m.values[f"high_exp_decay{s}"] = 1.0
        print(f"  High tail{s} REJECTED")
    return best_nll


def _test_inter_peak_group(m, gaus_amp1, gaus_amp2, best_nll):
    """Test peak1 high tail and peak2 low-side group together.

    These components overlap in the inter-peak region, so they must be tested
    jointly to avoid one absorbing the contribution of the other.
    Returns updated best_nll.
    """
    print("Testing inter-peak group (peak1 high tail + peak2 low-side)...")

    # Release peak1 high tail
    m.fixed["high_exp_frac1"] = False
    m.fixed["high_exp_decay1"] = False
    m.limits["high_exp_frac1"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["high_exp_decay1"] = (0.1, 50)
    m.values["high_exp_frac1"] = 0.15
    m.values["high_exp_decay1"] = 1.0

    # Release peak2 low-side group
    m.fixed["step_frac2"] = False
    m.limits["step_frac2"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.values["step_frac2"] = 0.1

    m.fixed["low_exp_frac2"] = False
    m.fixed["low_exp_decay2"] = False
    m.limits["low_exp_frac2"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["low_exp_decay2"] = (0.1, 50)
    m.values["low_exp_frac2"] = 0.15
    m.values["low_exp_decay2"] = 1.0

    m.fixed["low_lin_frac2"] = False
    m.fixed["low_lin_slope2"] = False
    m.limits["low_lin_frac2"] = (0, _MAX_SUBSIDIARY_FRAC)
    m.limits["low_lin_slope2"] = (-1, 1)
    m.values["low_lin_frac2"] = 0.15
    m.values["low_lin_slope2"] = 0.0

    m.migrad()
    print(
        f"  Inter-peak group NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
    )

    if m.fval < best_nll - 1:
        print("  Inter-peak group ACCEPTED")
        best_nll = m.fval
    else:
        print("  Inter-peak group REJECTED — re-fixing all")
        m.fixed["high_exp_frac1"] = True
        m.fixed["high_exp_decay1"] = True
        m.values["high_exp_frac1"] = 0.0
        m.values["high_exp_decay1"] = 1.0
        m.fixed["step_frac2"] = True
        m.values["step_frac2"] = 0.0
        m.fixed["low_exp_frac2"] = True
        m.fixed["low_exp_decay2"] = True
        m.values["low_exp_frac2"] = 0.0
        m.values["low_exp_decay2"] = 1.0
        m.fixed["low_lin_frac2"] = True
        m.fixed["low_lin_slope2"] = True
        m.values["low_lin_frac2"] = 0.0
        m.values["low_lin_slope2"] = 0.0

    return best_nll


def fit_double_peak(df_column,
                    fit_range_low,
                    fit_range_high,
                    expected_mu1,
                    expected_mu2,
                    cache_path=None):
    """Fit a double-peak model with automated component testing and pruning.

    Returns a Minuit object. If cache_path is given, results are cached.
    """
    if expected_mu1 > expected_mu2:
        expected_mu1, expected_mu2 = expected_mu2, expected_mu1

    data = np.asarray(df_column, dtype=np.float64)
    data = data[(data >= fit_range_low) & (data <= fit_range_high)]
    pdf = double_peak_pdf(fit_range_low, fit_range_high)
    c = cost.UnbinnedNLL(data, pdf)
    values, limits, fixed = estimate_double_peak_params(data,
                                                        fit_range_low,
                                                        fit_range_high,
                                                        mu1=expected_mu1,
                                                        mu2=expected_mu2)
    m = Minuit(c, **values)
    for k, v in limits.items():
        m.limits[k] = v
    for k, v in fixed.items():
        m.fixed[k] = v

    if cache_path and Path(cache_path).exists():
        _restore_minuit_state(m, cache_path)
        return m

    m.migrad()
    best_nll = m.fval
    gaus_amp1 = m.values["gaus_amp1"]
    gaus_amp2 = m.values["gaus_amp2"]
    print(f"Initial fit: NLL = {best_nll:.2f}")

    # Outer components (no inter-peak overlap)
    best_nll = _test_low_side_group(m, "1", gaus_amp1, best_nll)
    best_nll = _test_high_tail(m, "2", gaus_amp2, best_nll)
    # Inter-peak components tested jointly: peak1 high tail + peak2 low-side
    best_nll = _test_inter_peak_group(m, gaus_amp1, gaus_amp2, best_nll)

    print("Final fit with selected components...")
    m.migrad()
    m.hesse()

    # Ensure mu1 < mu2
    if m.values["mu1"] > m.values["mu2"]:
        print("Warning: mu1 > mu2 after fit, swapping peak parameters")
        for base in [
                "mu", "sigma", "gaus_amp", "step_frac", "low_exp_frac",
                "low_exp_decay", "low_lin_frac", "low_lin_slope",
                "high_exp_frac", "high_exp_decay"
        ]:
            v1, v2 = m.values[f"{base}1"], m.values[f"{base}2"]
            m.values[f"{base}1"] = v2
            m.values[f"{base}2"] = v1
            e1, e2 = m.errors[f"{base}1"], m.errors[f"{base}2"]
            m.errors[f"{base}1"] = e2
            m.errors[f"{base}2"] = e1
            f1, f2 = m.fixed[f"{base}1"], m.fixed[f"{base}2"]
            m.fixed[f"{base}1"] = f2
            m.fixed[f"{base}2"] = f1

    print(f"Final NLL = {m.fval:.2f}")
    print(m)
    if cache_path:
        _save_minuit_state(m, cache_path)
    return m
