"""Bin-integrated parametric model builders for simultaneous fits.

Each builder returns a callable ``f(edges, *params)`` that evaluates the model
integrated exactly over each bin (not evaluated at bin centers), consistent
with ``iminuit.cost.ExtendedBinnedNLL`` when ``extended=True``.

Parameter sharing between channels is by name: the builder's keyword
arguments set the parameter names exposed on the returned callable.

Returned callables carry their parameter names via ``__signature__``, which is
what :class:`iminuit.cost.ExtendedBinnedNLL` inspects.
"""

from __future__ import annotations

import inspect
from typing import Callable, Dict, Sequence, Tuple

import numpy as np
from scipy.special import erf, erfc, expit


def _with_param_names(core: Callable, *names: str) -> Callable:
    """Wrap ``core(edges, *args)`` so its signature exposes named params."""
    params: list[inspect.Parameter] = [
        inspect.Parameter("edges", inspect.Parameter.POSITIONAL_OR_KEYWORD),
    ]
    for name in names:
        params.append(
            inspect.Parameter(name, inspect.Parameter.POSITIONAL_OR_KEYWORD))
    signature = inspect.Signature(params)

    def wrapper(edges: np.ndarray, *values: float) -> np.ndarray:
        return core(edges, *values)

    wrapper.__signature__ = signature
    wrapper.__name__ = getattr(core, "__name__", "model")
    return wrapper


def gaussian_peak(
    mu: str = "mu",
    sigma: str = "sigma",
    nsig: str = "nsig",
) -> Callable:
    """Gaussian peak, bin-integrated via ``erf``.

    Parameters
    ----------
    mu, sigma, nsig : str
        Names of the peak center, width, and yield parameters on the
        returned callable. Rename to share or isolate across channels.

    Returns
    -------
    callable
        ``f(edges, mu, sigma, nsig) -> ndarray`` of length ``len(edges)-1``.
        ``nsig`` is the yield integrated over :math:`(-\\infty, +\\infty)`.
    """
    def core(edges: np.ndarray, mu_val: float, sigma_val: float,
             nsig_val: float) -> np.ndarray:
        z = (np.asarray(edges) - mu_val) / (np.sqrt(2.0) * sigma_val)
        cdf = 0.5 * (1.0 + erf(z))
        return nsig_val * np.diff(cdf)

    return _with_param_names(core, mu, sigma, nsig)


def exponential_bkg(
    reference_range: Tuple[float, float],
    tau: str = "tau",
    nbkg: str = "nbkg",
) -> Callable:
    """Exponential background :math:`\\exp(-x/\\tau)`, bin-integrated.

    Parameters
    ----------
    reference_range : (float, float)
        ``(e_lo, e_hi)`` range over which ``nbkg`` is defined as the integral.
    tau : str
        Name of the decay constant parameter.
    nbkg : str
        Name of the yield parameter (integral over ``reference_range``).

    Returns
    -------
    callable
        ``f(edges, tau, nbkg) -> ndarray`` of per-bin expected counts.
    """
    e_lo, e_hi = float(reference_range[0]), float(reference_range[1])

    def core(edges: np.ndarray, tau_val: float,
             nbkg_val: float) -> np.ndarray:
        edges_arr = np.asarray(edges, dtype=np.float64)
        shifted_hi = e_hi - e_lo
        shifted_edges = edges_arr - e_lo
        ref_integral = 1.0 - np.exp(-shifted_hi / tau_val)
        bin_integrals = (np.exp(-shifted_edges[:-1] / tau_val)
                         - np.exp(-shifted_edges[1:] / tau_val))
        return nbkg_val * bin_integrals / ref_integral

    return _with_param_names(core, tau, nbkg)


def linear_bkg(
    reference_range: Tuple[float, float],
    slope: str = "slope",
    nbkg: str = "nbkg",
) -> Callable:
    """Linear background :math:`1 + s \\cdot (x - x_\\mathrm{mid})`, bin-integrated.

    Parametrized about the midpoint of ``reference_range`` so that the
    reference-range integral is independent of slope at nominal; ``nbkg`` is
    therefore the integral over ``reference_range`` regardless of slope.

    Parameters
    ----------
    reference_range : (float, float)
    slope : str
    nbkg : str

    Returns
    -------
    callable
        ``f(edges, slope, nbkg) -> ndarray``.
    """
    e_lo, e_hi = float(reference_range[0]), float(reference_range[1])
    e_mid = 0.5 * (e_lo + e_hi)
    ref_width = e_hi - e_lo

    def core(edges: np.ndarray, slope_val: float,
             nbkg_val: float) -> np.ndarray:
        edges_arr = np.asarray(edges, dtype=np.float64)
        dx = edges_arr[1:] - edges_arr[:-1]
        mid = 0.5 * (edges_arr[1:] + edges_arr[:-1])
        density = 1.0 + slope_val * (mid - e_mid)
        return (nbkg_val / ref_width) * density * dx

    return _with_param_names(core, slope, nbkg)


def polynomial_bkg(
    reference_range: Tuple[float, float],
    coeff_names: Sequence[str],
    nbkg: str = "nbkg",
) -> Callable:
    """Polynomial background :math:`1 + \\sum_{i\\ge 1} c_i (x - x_\\mathrm{mid})^i`.

    Analytically integrated per bin; normalized so the integral over
    ``reference_range`` equals ``nbkg`` regardless of coefficient values.

    Parameters
    ----------
    reference_range : (float, float)
    coeff_names : sequence of str
        Names for coefficients :math:`c_1, c_2, \\ldots`. The constant
        (intercept) is fixed at 1 and absorbed into ``nbkg``; pass an empty
        sequence for a flat background.
    nbkg : str

    Returns
    -------
    callable
        ``f(edges, *coeffs, nbkg) -> ndarray``.
    """
    e_lo, e_hi = float(reference_range[0]), float(reference_range[1])
    e_mid = 0.5 * (e_lo + e_hi)
    n_coeffs = len(coeff_names)

    def antideriv(x: np.ndarray, coeffs: Sequence[float]) -> np.ndarray:
        x_shift = x - e_mid
        total = x_shift.copy() if isinstance(x_shift, np.ndarray) else x_shift
        for i in range(n_coeffs):
            power = i + 2
            total = total + coeffs[i] * x_shift**power / power
        return total

    def core(edges: np.ndarray, *args: float) -> np.ndarray:
        coeffs = args[:n_coeffs]
        nbkg_val = args[n_coeffs]
        edges_arr = np.asarray(edges, dtype=np.float64)
        ref_integral = (antideriv(np.array(e_hi), coeffs)
                        - antideriv(np.array(e_lo), coeffs))
        bin_integrals = (antideriv(edges_arr[1:], coeffs)
                         - antideriv(edges_arr[:-1], coeffs))
        return nbkg_val * bin_integrals / ref_integral

    return _with_param_names(core, *coeff_names, nbkg)


def step_shelf(
    mu: str = "mu",
    sigma: str = "sigma",
    amplitude: str = "A_step",
) -> Callable:
    """Fermi-Dirac-squared step shelf, bin-integrated.

    Density form (matches :cpp:class:`FittingUtils`):

    .. math::

       f(x) = \\frac{A}{\\left(1 + e^{(x - \\mu)/\\sigma}\\right)^2}

    At :math:`x \\to -\\infty` the density tends to ``A``; at :math:`x = \\mu`
    the density is ``A/4``; at :math:`x \\to +\\infty` the density tends to 0.
    Used to model events in which part of the photon energy escaped the
    active volume, smeared by the detector resolution.

    The antiderivative with respect to ``z = (x - mu) / sigma`` is
    :math:`F(z) = -\\log(1 + e^{-z}) + 1 / (1 + e^z)`; bin integrals are
    ``A * sigma * (F(z_hi) - F(z_lo))``.

    Parameters
    ----------
    mu, sigma, amplitude : str
        Parameter names on the returned callable. Use the same ``mu``/
        ``sigma`` as a sibling gaussian to tie the shelf to that peak.

    Returns
    -------
    callable
        ``f(edges, mu, sigma, amplitude) -> ndarray``.
    """
    def core(edges: np.ndarray, mu_val: float, sigma_val: float,
             amplitude_val: float) -> np.ndarray:
        edges_arr = np.asarray(edges, dtype=np.float64)
        z = (edges_arr - mu_val) / sigma_val
        antideriv = -np.logaddexp(0.0, -z) + expit(-z)
        return amplitude_val * sigma_val * np.diff(antideriv)

    return _with_param_names(core, mu, sigma, amplitude)


def _exp_tail_antideriv(
    y: np.ndarray,
    sigma_val: float,
    tau_val: float,
) -> np.ndarray:
    """Antiderivative of ``exp(y/tau) * erfc(y/(sigma*sqrt(2)))`` w.r.t. ``y``."""
    sqrt2_sigma = np.sqrt(2.0) * sigma_val
    offset = sigma_val * sigma_val / tau_val
    gauss_corr = np.exp(sigma_val * sigma_val / (2.0 * tau_val * tau_val))
    return tau_val * (np.exp(y / tau_val) * erfc(y / sqrt2_sigma)
                      + gauss_corr * erf((y - offset) / sqrt2_sigma))


def low_exp_tail(
    mu: str = "mu",
    sigma: str = "sigma",
    amplitude: str = "A_low_exp",
    decay: str = "tau_low_exp",
) -> Callable:
    """Low-energy exponential tail, bin-integrated.

    Density form (matches :cpp:class:`FittingUtils`):

    .. math::

       f(x) = A \\, e^{(x - \\mu)/\\tau} \\,
              \\operatorname{erfc}\\!\\left(\\frac{x - \\mu}{\\sigma \\sqrt 2}\\right)

    Rises exponentially below the peak and is suppressed by ``erfc`` above
    it (smeared step at the peak position). Models incomplete charge
    collection / trapping.

    Parameters
    ----------
    mu, sigma, amplitude, decay : str

    Returns
    -------
    callable
        ``f(edges, mu, sigma, amplitude, decay) -> ndarray``.
    """
    def core(edges: np.ndarray, mu_val: float, sigma_val: float,
             amplitude_val: float, tau_val: float) -> np.ndarray:
        y = np.asarray(edges, dtype=np.float64) - mu_val
        return amplitude_val * np.diff(
            _exp_tail_antideriv(y, sigma_val, tau_val))

    return _with_param_names(core, mu, sigma, amplitude, decay)


def high_exp_tail(
    mu: str = "mu",
    sigma: str = "sigma",
    amplitude: str = "A_high_exp",
    decay: str = "tau_high_exp",
) -> Callable:
    """High-energy exponential tail, bin-integrated.

    Mirror of :func:`low_exp_tail` about ``mu``:

    .. math::

       f(x) = A \\, e^{(\\mu - x)/\\tau} \\,
              \\operatorname{erfc}\\!\\left(\\frac{\\mu - x}{\\sigma \\sqrt 2}\\right)

    Used to capture pileup-induced tails above the photopeak.

    Parameters
    ----------
    mu, sigma, amplitude, decay : str

    Returns
    -------
    callable
        ``f(edges, mu, sigma, amplitude, decay) -> ndarray``.
    """
    def core(edges: np.ndarray, mu_val: float, sigma_val: float,
             amplitude_val: float, tau_val: float) -> np.ndarray:
        z = mu_val - np.asarray(edges, dtype=np.float64)
        return -amplitude_val * np.diff(
            _exp_tail_antideriv(z, sigma_val, tau_val))

    return _with_param_names(core, mu, sigma, amplitude, decay)


def low_lin_tail(
    mu: str = "mu",
    sigma: str = "sigma",
    amplitude: str = "A_low_lin",
    slope: str = "s_low_lin",
) -> Callable:
    """Low-energy linear tail, bin-integrated.

    Density form (matches :cpp:class:`FittingUtils`):

    .. math::

       f(x) = A \\, \\left[1 + s (x - \\mu)\\right] \\,
              \\operatorname{erfc}\\!\\left(\\frac{x - \\mu}{\\sigma \\sqrt 2}\\right)

    Captures asymmetric tailing below the peak not well described by a
    single exponential.

    Parameters
    ----------
    mu, sigma, amplitude, slope : str

    Returns
    -------
    callable
        ``f(edges, mu, sigma, amplitude, slope) -> ndarray``.
    """
    def core(edges: np.ndarray, mu_val: float, sigma_val: float,
             amplitude_val: float, slope_val: float) -> np.ndarray:
        y = np.asarray(edges, dtype=np.float64) - mu_val
        antideriv = _low_lin_antideriv(y, sigma_val, slope_val)
        return amplitude_val * np.diff(antideriv)

    return _with_param_names(core, mu, sigma, amplitude, slope)


def _low_lin_antideriv(
    y: np.ndarray,
    sigma_val: float,
    slope_val: float,
) -> np.ndarray:
    """Antiderivative of ``(1 + s*y) * erfc(y/(sigma*sqrt(2)))`` w.r.t. ``y``."""
    beta = 1.0 / (np.sqrt(2.0) * sigma_val)
    by = beta * y
    erfc_by = erfc(by)
    erf_by = erf(by)
    gauss = np.exp(-(by * by))
    flat = y * erfc_by - gauss / (beta * np.sqrt(np.pi))
    lin = (0.5 * y * y * erfc_by
           - y * gauss / (2.0 * beta * np.sqrt(np.pi))
           + erf_by / (4.0 * beta * beta))
    return flat + slope_val * lin


def photopeak(
    prefix: str,
    mu: str = "mu",
    sigma: str = "sigma",
    nsig: str = "nsig",
    low_exp: bool = False,
    low_lin: bool = False,
    high_exp: bool = False,
    step: bool = False,
) -> Callable:
    """Build a gaussian peak with optional shelf/tail components.

    Convenience wrapper around :func:`compose` for the FittingUtils-style
    photopeak: gaussian plus any subset of {step shelf, low-energy
    exponential tail, low-energy linear tail, high-energy exponential
    tail}. All tail components share the gaussian's ``mu`` and ``sigma``
    by name; tail-specific parameters (amplitude, decay, slope) are
    prefixed with ``prefix`` so that multiple photopeaks in the same fit
    keep distinct parameter names.

    Parameters
    ----------
    prefix : str
        Prepended to tail-specific parameter names. E.g. ``prefix="pb1_"``
        produces ``pb1_A_low_exp``, ``pb1_tau_low_exp``, etc.
    mu, sigma, nsig : str
        Parameter names for the central gaussian (not prefixed).
    low_exp, low_lin, high_exp, step : bool
        Enable each optional component.

    Returns
    -------
    callable
        A composed callable with a ``.components`` dict whose keys are
        ``"gauss"`` plus any of ``"step"``, ``"low_exp"``, ``"low_lin"``,
        ``"high_exp"``.
    """
    components: Dict[str, Callable] = {
        "gauss": gaussian_peak(mu=mu, sigma=sigma, nsig=nsig),
    }
    if step:
        components["step"] = step_shelf(
            mu=mu, sigma=sigma, amplitude=f"{prefix}A_step")
    if low_exp:
        components["low_exp"] = low_exp_tail(
            mu=mu,
            sigma=sigma,
            amplitude=f"{prefix}A_low_exp",
            decay=f"{prefix}tau_low_exp",
        )
    if low_lin:
        components["low_lin"] = low_lin_tail(
            mu=mu,
            sigma=sigma,
            amplitude=f"{prefix}A_low_lin",
            slope=f"{prefix}s_low_lin",
        )
    if high_exp:
        components["high_exp"] = high_exp_tail(
            mu=mu,
            sigma=sigma,
            amplitude=f"{prefix}A_high_exp",
            decay=f"{prefix}tau_high_exp",
        )
    return compose(**components)


def compose(**components: Callable) -> Callable:
    """Sum named component callables into a total model with overlay metadata.

    All components must share the ``f(edges, *params)`` convention. Parameter
    sharing between components is by name: a parameter appearing in multiple
    component signatures is a single fit parameter.

    Parameters
    ----------
    **components
        Keyword-named component callables. The keys become the labels used by
        :meth:`SimultaneousFit.plot` for overlay graphs.

    Returns
    -------
    callable
        A model callable with attributes:
          * ``.components`` : ``dict[str, callable]`` of the input components.
          * ``__signature__`` : union of component parameter names in insertion order.
    """
    if not components:
        raise ValueError("compose() requires at least one component")

    all_params: list[str] = []
    seen: set[str] = set()
    component_params: Dict[str, list[str]] = {}
    for comp_name, comp in components.items():
        sig = inspect.signature(comp)
        pnames = [p.name for p in list(sig.parameters.values())[1:]]
        component_params[comp_name] = pnames
        for p in pnames:
            if p not in seen:
                all_params.append(p)
                seen.add(p)

    def core(edges: np.ndarray, *values: float) -> np.ndarray:
        arg_map = {name: val for name, val in zip(all_params, values)}
        result: np.ndarray | None = None
        for comp_name, comp in components.items():
            comp_args = [arg_map[p] for p in component_params[comp_name]]
            contribution = comp(edges, *comp_args)
            result = contribution if result is None else result + contribution
        return result

    total = _with_param_names(core, *all_params)
    total.components = dict(components)
    return total

