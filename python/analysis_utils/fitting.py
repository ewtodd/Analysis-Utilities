import numpy as np
from scipy.special import erfc

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
    z = (x - mu) / sigma
    denom = (1.0 + np.exp(z))**2
    return np.where(denom < 1e-100, 0.0, step_amplitude / denom)


def _low_tail(x, mu, sigma, exp_amp, exp_decay, lin_amp, lin_slope):
    """Low-energy tail (exponential + linear) * erfc.

    Mirrors FittingFunctions::LowTail.
    """
    if sigma <= 0.0 or (exp_amp == 0.0 and lin_amp == 0.0):
        return np.zeros_like(x, dtype=np.float64)

    y = x - mu
    exp_term = exp_amp * np.exp(y / exp_decay) if exp_amp != 0.0 else 0.0
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
    return exp_amp * np.exp(y / exp_decay) * erfc(y / (_SQRT2 * sigma))


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
