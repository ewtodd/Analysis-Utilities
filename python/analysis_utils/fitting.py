"""PDF component functions for unbinned likelihood fitting.

Ports the C++ FittingFunctions namespace from FittingUtils.cpp.
All functions accept numpy arrays or scalars for x; parameters are scalars.
"""

from functools import partial

import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
from iminuit import Minuit, cost

_SQRT2 = np.sqrt(2.0)


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

    # Step
    if use_step:
        values[f"step_amp{suffix}"] = peak_height * 0.1
        limits[f"step_amp{suffix}"] = (0, peak_height)
    else:
        values[f"step_amp{suffix}"] = 0.0
        fixed[f"step_amp{suffix}"] = True

    # Low exponential tail
    if use_low_exp_tail:
        values[f"low_exp_amp{suffix}"] = min(peak_height * 0.15,
                                             peak_height * 0.25)
        limits[f"low_exp_amp{suffix}"] = (0, peak_height)
        values[f"low_exp_decay{suffix}"] = 1.0
        limits[f"low_exp_decay{suffix}"] = (0.1, 50)
    else:
        values[f"low_exp_amp{suffix}"] = 0.0
        values[f"low_exp_decay{suffix}"] = 1.0
        fixed[f"low_exp_amp{suffix}"] = True
        fixed[f"low_exp_decay{suffix}"] = True

    # Low linear tail
    if use_low_lin_tail:
        values[f"low_lin_amp{suffix}"] = min(peak_height * 0.15,
                                             peak_height * 0.25)
        limits[f"low_lin_amp{suffix}"] = (0, peak_height)
        values[f"low_lin_slope{suffix}"] = 0.0
        limits[f"low_lin_slope{suffix}"] = (-0.5 * bkg_estimate / range_width,
                                            0.5 * bkg_estimate / range_width)
    else:
        values[f"low_lin_amp{suffix}"] = 0.0
        values[f"low_lin_slope{suffix}"] = 0.0
        fixed[f"low_lin_amp{suffix}"] = True
        fixed[f"low_lin_slope{suffix}"] = True

    # High exponential tail
    if use_high_exp_tail:
        values[f"high_exp_amp{suffix}"] = min(peak_height * 0.15,
                                              peak_height * 0.25)
        limits[f"high_exp_amp{suffix}"] = (0, peak_height)
        values[f"high_exp_decay{suffix}"] = 1.0
        limits[f"high_exp_decay{suffix}"] = (0.1, 50)
    else:
        values[f"high_exp_amp{suffix}"] = 0.0
        values[f"high_exp_decay{suffix}"] = 1.0
        fixed[f"high_exp_amp{suffix}"] = True
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


def _single_peak_pdf(fit_range_low, fit_range_high, quad_limit, x, mu, sigma,
                     gaus_amp, step_amp, low_exp_amp, low_exp_decay,
                     low_lin_amp, low_lin_slope, high_exp_amp, high_exp_decay,
                     bkg_constant, lin_bkg_slope):
    args = (mu, sigma, gaus_amp, step_amp, low_exp_amp, low_exp_decay,
            low_lin_amp, low_lin_slope, high_exp_amp, high_exp_decay,
            bkg_constant, lin_bkg_slope)
    norm, _ = quad(_peak_function,
                   fit_range_low,
                   fit_range_high,
                   args=args,
                   limit=quad_limit,
                   epsabs=1e-10,
                   epsrel=1e-8)
    return _peak_function(x, *args) / norm


def _double_peak_pdf(fit_range_low, fit_range_high, quad_limit, x, mu1, sigma1,
                     gaus_amp1, step_amp1, low_exp_amp1, low_exp_decay1,
                     low_lin_amp1, low_lin_slope1, high_exp_amp1,
                     high_exp_decay1, mu2, sigma2, gaus_amp2, step_amp2,
                     low_exp_amp2, low_exp_decay2, low_lin_amp2,
                     low_lin_slope2, high_exp_amp2, high_exp_decay2,
                     bkg_constant, lin_bkg_slope):
    args = (mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1, low_exp_decay1,
            low_lin_amp1, low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
            sigma2, gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
            low_lin_amp2, low_lin_slope2, high_exp_amp2, high_exp_decay2,
            bkg_constant, lin_bkg_slope)
    norm, _ = quad(_double_peak_function,
                   fit_range_low,
                   fit_range_high,
                   args=args,
                   limit=quad_limit,
                   epsabs=1e-10,
                   epsrel=1e-8)
    return _double_peak_function(x, *args) / norm


def _triple_peak_pdf(
        fit_range_low, fit_range_high, quad_limit, x, mu1, sigma1, gaus_amp1,
        step_amp1, low_exp_amp1, low_exp_decay1, low_lin_amp1, low_lin_slope1,
        high_exp_amp1, high_exp_decay1, mu2, sigma2, gaus_amp2, step_amp2,
        low_exp_amp2, low_exp_decay2, low_lin_amp2, low_lin_slope2,
        high_exp_amp2, high_exp_decay2, mu3, sigma3, gaus_amp3, step_amp3,
        low_exp_amp3, low_exp_decay3, low_lin_amp3, low_lin_slope3,
        high_exp_amp3, high_exp_decay3, bkg_constant, lin_bkg_slope):
    args = (mu1, sigma1, gaus_amp1, step_amp1, low_exp_amp1, low_exp_decay1,
            low_lin_amp1, low_lin_slope1, high_exp_amp1, high_exp_decay1, mu2,
            sigma2, gaus_amp2, step_amp2, low_exp_amp2, low_exp_decay2,
            low_lin_amp2, low_lin_slope2, high_exp_amp2, high_exp_decay2, mu3,
            sigma3, gaus_amp3, step_amp3, low_exp_amp3, low_exp_decay3,
            low_lin_amp3, low_lin_slope3, high_exp_amp3, high_exp_decay3,
            bkg_constant, lin_bkg_slope)
    norm, _ = quad(_triple_peak_function,
                   fit_range_low,
                   fit_range_high,
                   args=args,
                   limit=quad_limit,
                   epsabs=1e-10,
                   epsrel=1e-8)
    return _triple_peak_function(x, *args) / norm


def single_peak_pdf(fit_range_low, fit_range_high, quad_limit=50):
    """Return a normalized single-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_single_peak_pdf, fit_range_low, fit_range_high, quad_limit)


def double_peak_pdf(fit_range_low, fit_range_high, quad_limit=50):
    """Return a normalized double-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_double_peak_pdf, fit_range_low, fit_range_high, quad_limit)


def triple_peak_pdf(fit_range_low, fit_range_high, quad_limit=50):
    """Return a normalized triple-peak PDF for use with cost.UnbinnedNLL."""
    return partial(_triple_peak_pdf, fit_range_low, fit_range_high, quad_limit)


def fit_single_peak(df_column, fit_range_low, fit_range_high, expected_mu):
    """
    Return a Minuit object after fitting with UnbinnedNLL cost function is completed. 
    """
    pdf = single_peak_pdf(fit_range_low, fit_range_high)
    c = cost.UnbinnedNLL(df_column, pdf)
    values, limits, fixed = estimate_peak_params(df_column,
                                                 fit_range_low,
                                                 fit_range_high,
                                                 mu=expected_mu)
    m = Minuit(c, **values)
    for k, v in limits.items():
        m.limits[k] = v
    for k, v in fixed.items():
        m.fixed[k] = v

    m.migrad()
    best_nll = m.fval
    gaus_amp = m.values["gaus_amp"]
    print(f"Initial fit: NLL = {best_nll:.2f}")

    print("Testing low-side group (step + low exp tail + low lin tail)...")
    m.fixed["step_amp"] = False
    m.limits["step_amp"] = (0, gaus_amp * 2)
    m.values["step_amp"] = gaus_amp * 0.1

    m.fixed["low_exp_amp"] = False
    m.fixed["low_exp_decay"] = False
    m.limits["low_exp_amp"] = (0, gaus_amp * 2)
    m.limits["low_exp_decay"] = (0.1, 50)
    m.values["low_exp_amp"] = gaus_amp * 0.15
    m.values["low_exp_decay"] = 1.0

    m.fixed["low_lin_amp"] = False
    m.fixed["low_lin_slope"] = False
    m.limits["low_lin_amp"] = (0, gaus_amp * 2)
    m.limits["low_lin_slope"] = (-1, 1)
    m.values["low_lin_amp"] = gaus_amp * 0.15
    m.values["low_lin_slope"] = 0.0

    m.migrad()
    print(
        f"  Low-side group NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
    )

    if m.fval < best_nll - 1:
        print("  Low-side group ACCEPTED — pruning individual components...")
        best_nll = m.fval

        saved = m.values["step_amp"]
        m.fixed["step_amp"] = True
        m.values["step_amp"] = 0.0
        m.migrad()
        print(
            f"  Without step: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Step PRUNED")
        else:
            m.fixed["step_amp"] = False
            m.values["step_amp"] = saved
            m.migrad()
            print("  Step KEPT")

        saved_amp = m.values["low_exp_amp"]
        saved_decay = m.values["low_exp_decay"]
        m.fixed["low_exp_amp"] = True
        m.fixed["low_exp_decay"] = True
        m.values["low_exp_amp"] = 0.0
        m.values["low_exp_decay"] = 1.0
        m.migrad()
        print(
            f"  Without low exp tail: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Low exp tail PRUNED")
        else:
            m.fixed["low_exp_amp"] = False
            m.fixed["low_exp_decay"] = False
            m.values["low_exp_amp"] = saved_amp
            m.values["low_exp_decay"] = saved_decay
            m.migrad()
            print("  Low exp tail KEPT")

        saved_amp = m.values["low_lin_amp"]
        saved_slope = m.values["low_lin_slope"]
        m.fixed["low_lin_amp"] = True
        m.fixed["low_lin_slope"] = True
        m.values["low_lin_amp"] = 0.0
        m.values["low_lin_slope"] = 0.0
        m.migrad()
        print(
            f"  Without low lin tail: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Low lin tail PRUNED")
        else:
            m.fixed["low_lin_amp"] = False
            m.fixed["low_lin_slope"] = False
            m.values["low_lin_amp"] = saved_amp
            m.values["low_lin_slope"] = saved_slope
            m.migrad()
            print("  Low lin tail KEPT")
    else:
        print("  Low-side group REJECTED — re-fixing all")
        m.fixed["step_amp"] = True
        m.values["step_amp"] = 0.0
        m.fixed["low_exp_amp"] = True
        m.fixed["low_exp_decay"] = True
        m.values["low_exp_amp"] = 0.0
        m.values["low_exp_decay"] = 1.0
        m.fixed["low_lin_amp"] = True
        m.fixed["low_lin_slope"] = True
        m.values["low_lin_amp"] = 0.0
        m.values["low_lin_slope"] = 0.0

    print("Testing high exponential tail...")
    m.fixed["high_exp_amp"] = False
    m.fixed["high_exp_decay"] = False
    m.limits["high_exp_amp"] = (0, gaus_amp * 2)
    m.limits["high_exp_decay"] = (0.1, 50)
    m.values["high_exp_amp"] = gaus_amp * 0.15
    m.values["high_exp_decay"] = 1.0
    m.migrad()
    print(f"  High tail NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})")
    if m.fval < best_nll - 1:
        best_nll = m.fval
        print("  High tail ACCEPTED")
    else:
        m.fixed["high_exp_amp"] = True
        m.fixed["high_exp_decay"] = True
        m.values["high_exp_amp"] = 0.0
        m.values["high_exp_decay"] = 1.0
        print("  High tail REJECTED")

    print(f"Final fit with selected components...")
    m.migrad()
    m.hesse()
    print(f"Final NLL = {m.fval:.2f}")
    print(m)
    return m
