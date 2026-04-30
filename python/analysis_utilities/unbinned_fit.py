"""Unbinned extended maximum-likelihood fit on event-level energy data.

Parallel to :class:`SimultaneousFit` for the binned case, but a single
channel of 1D events. Built on :class:`iminuit.cost.ExtendedUnbinnedNLL`:
the model callable returns ``(total_yield, density)`` over a fit window
``(x_lo, x_hi)``, where ``density`` is normalized to integrate to 1 over
that window.

Use this as a cross-check / narrow-range alternative to binned simultaneous
fits — e.g. fit only the Ge peak with an exponential continuum, sweep the
upper fit-range cutoff, and quote the spread as a fit-range systematic.
"""

from __future__ import annotations

import array as _array
import inspect
from dataclasses import dataclass
from typing import (Any, Callable, Dict, Iterable, List, Optional, Sequence,
                    Tuple)

import numpy as np
import pandas as pd


@dataclass
class _ParameterSpec:
    value: Optional[float] = None
    limits: Optional[Tuple[Optional[float], Optional[float]]] = None
    fixed: bool = False
    error: Optional[float] = None


@dataclass
class _ConstraintSpec:
    name: str
    mean: float
    sigma: float


@dataclass(frozen=True)
class UnbinnedFitResult:
    """Structured unbinned-fit output.

    Mirrors :class:`analysis_utilities.simultaneous_fit.FitResult` but for
    the single-channel unbinned case. See that class for field semantics.
    """

    valid: bool
    values: Dict[str, float]
    errors: Dict[str, float]
    minos: Dict[str, Tuple[float, float]]
    corr: pd.DataFrame
    nll: float
    chi2: float
    ndof: int
    x_range: Tuple[float, float]
    at_limit: List[str]
    minuit: Any


class UnbinnedFit:
    """Event-level extended maximum-likelihood fit on a fixed energy window.

    Parameters
    ----------
    name : str, optional
        Output-name prefix for plots.

    Notes
    -----
    Call :meth:`set_events`, register parameters, then :meth:`fit`. The
    model callable supplied to :meth:`set_events` must have signature
    ``model(x, *params) -> (total_yield, density)`` where ``density``
    integrates to 1 over the fit range — the
    :func:`analysis_utilities.fit_models.gaussian_peak_density` and
    :func:`analysis_utilities.fit_models.exponential_bkg_density`
    helpers produce components that compose additively.
    """

    def __init__(self, name: str = "unbinned") -> None:
        self.name: str = name
        self._events: Optional[np.ndarray] = None
        self._x_range: Optional[Tuple[float, float]] = None
        self._model: Optional[Callable] = None
        self._components: Optional[Dict[str, Callable]] = None
        self._background_component: Optional[str] = None
        self._param_specs: Dict[str, _ParameterSpec] = {}
        self._constraints: List[_ConstraintSpec] = []

    def set_events(
        self,
        events: np.ndarray,
        x_range: Tuple[float, float],
        model: Callable,
        components: Optional[Dict[str, Callable]] = None,
        background_component: Optional[str] = None,
    ) -> None:
        """Register the 1D event energies and model.

        Parameters
        ----------
        events : array-like
            1D array of event energies. Events outside ``x_range`` are
            dropped on registration.
        x_range : tuple(float, float)
            Fit window ``(x_lo, x_hi)``. Density must integrate to 1 over
            this window; total_yield is the expected event count inside
            the window.
        model : callable
            ``model(x, *params) -> (total_yield, density)`` — signature
            expected by :class:`iminuit.cost.ExtendedUnbinnedNLL`. Param
            names taken from ``inspect.signature``.
        components : dict[str, callable], optional
            Named sub-model callables for plot overlays. Each
            ``sub(x, *params) -> (yield, density)``.
        background_component : str, optional
            Key within ``components`` identifying the background; other
            components are drawn summed with it (FittingUtils convention).
        """
        x_lo, x_hi = float(x_range[0]), float(x_range[1])
        events_arr = np.asarray(events, dtype=np.float64)
        mask = (events_arr >= x_lo) & (events_arr <= x_hi)
        self._events = events_arr[mask].copy()
        self._x_range = (x_lo, x_hi)
        self._model = model
        if components is None:
            components = getattr(model, "components", None)
        self._components = components
        self._background_component = background_component

    def set_parameter(
        self,
        name: str,
        value: Optional[float] = None,
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        fixed: bool = False,
        error: Optional[float] = None,
    ) -> None:
        self._param_specs[name] = _ParameterSpec(
            value=value,
            limits=limits,
            fixed=fixed,
            error=error,
        )

    def add_constraint(self, name: str, mean: float, sigma: float) -> None:
        self._constraints.append(
            _ConstraintSpec(name=name, mean=float(mean), sigma=float(sigma)))

    def fit(
        self,
        minos: Optional[Iterable[str] | bool] = None,
        strategy: int = 1,
        tol: Optional[float] = None,
        print_level: int = 0,
        migrad_ncall: Optional[int] = None,
        migrad_iterate: int = 5,
    ) -> UnbinnedFitResult:
        """Run migrad + hesse (+ optional minos). See
        :meth:`SimultaneousFit.fit` for argument semantics."""
        from iminuit import Minuit
        from iminuit.cost import ExtendedUnbinnedNLL, NormalConstraint

        if self._events is None or self._model is None:
            raise RuntimeError("set_events() has not been called")

        cost: Any = ExtendedUnbinnedNLL(self._events, self._model)

        for con in self._constraints:
            nc = NormalConstraint(con.name, con.mean, con.sigma)
            cost = cost + nc

        all_params: List[str] = list(cost._parameters.keys())

        initial_values: Dict[str, float] = {}
        for pname in all_params:
            spec = self._param_specs.get(pname)
            if spec is None or spec.value is None:
                import warnings
                warnings.warn(
                    f"Parameter '{pname}' has no initial value; using 1.0",
                    stacklevel=2)
                initial_values[pname] = 1.0
            else:
                initial_values[pname] = float(spec.value)

        m = Minuit(cost, **initial_values)
        m.strategy = strategy
        m.print_level = print_level
        if tol is not None:
            m.tol = tol

        for pname in all_params:
            spec = self._param_specs.get(pname)
            if spec is None:
                continue
            if spec.limits is not None:
                m.limits[pname] = spec.limits
            if spec.fixed:
                m.fixed[pname] = True
            if spec.error is not None:
                m.errors[pname] = spec.error

        m.migrad(ncall=migrad_ncall, iterate=migrad_iterate)
        m.hesse()

        at_limit: List[str] = []
        for p in m.params:
            if p.is_fixed:
                continue
            scale = max(abs(p.value), 1.0)
            if (p.has_lower_limit
                    and abs(p.value - p.lower_limit) < 1e-4 * scale):
                at_limit.append(f"{p.name}@lower={p.lower_limit:g}")
            if (p.has_upper_limit
                    and abs(p.value - p.upper_limit) < 1e-4 * scale):
                at_limit.append(f"{p.name}@upper={p.upper_limit:g}")

        if not m.valid:
            import warnings
            warnings.warn(
                f"Minuit did not reach a valid minimum (fval={m.fval:.3f}, "
                f"edm={m.fmin.edm:.3g}). Parameters at limit: {at_limit}. "
                "Minos (if requested) will be skipped.",
                stacklevel=2)

        minos_results: Dict[str, Tuple[float, float]] = {}
        if minos is not None and minos is not False and m.valid:
            if minos is True:
                m.minos()
                for pname in all_params:
                    if not m.fixed[pname] and pname in m.merrors:
                        merr = m.merrors[pname]
                        minos_results[pname] = (float(merr.lower),
                                                 float(merr.upper))
            else:
                for pname in minos:
                    if pname not in all_params:
                        raise ValueError(
                            f"Minos requested for unknown parameter '{pname}'")
                    m.minos(pname)
                    merr = m.merrors[pname]
                    minos_results[pname] = (float(merr.lower),
                                             float(merr.upper))

        values: Dict[str, float] = {n: float(m.values[n]) for n in all_params}
        errors: Dict[str, float] = {n: float(m.errors[n]) for n in all_params}

        free_params: List[str] = [n for n in all_params if not m.fixed[n]]
        if m.covariance is not None and len(free_params) > 0:
            n_free = len(free_params)
            corr_arr = np.empty((n_free, n_free), dtype=np.float64)
            for i in range(n_free):
                for j in range(n_free):
                    ni, nj = free_params[i], free_params[j]
                    cov_ii = m.covariance[ni, ni]
                    cov_jj = m.covariance[nj, nj]
                    cov_ij = m.covariance[ni, nj]
                    if cov_ii > 0 and cov_jj > 0:
                        corr_arr[i, j] = cov_ij / np.sqrt(cov_ii * cov_jj)
                    else:
                        corr_arr[i, j] = np.nan
            corr = pd.DataFrame(corr_arr,
                                 index=free_params,
                                 columns=free_params)
        else:
            corr = pd.DataFrame()

        x_lo, x_hi = self._x_range
        auto_nbins = max(20, int(np.sqrt(len(self._events))))
        hist_edges = np.linspace(x_lo, x_hi, auto_nbins + 1)
        counts, _ = np.histogram(self._events, bins=hist_edges)

        param_names = _model_params(self._model)
        param_values = [values[p] for p in param_names]
        _, rate_at_edges = self._model(hist_edges, *param_values)
        bin_widths = np.diff(hist_edges)
        expected = 0.5 * (rate_at_edges[:-1] + rate_at_edges[1:]) * bin_widths

        with np.errstate(divide="ignore", invalid="ignore"):
            pulls = np.where(
                expected > 0,
                (counts - expected) / np.sqrt(expected),
                0.0,
            )
        chi2 = float(np.sum(pulls * pulls))
        ndof = int(len(counts) - len(free_params))

        return UnbinnedFitResult(
            valid=bool(m.valid and m.covariance is not None),
            values=values,
            errors=errors,
            minos=minos_results,
            corr=corr,
            nll=float(m.fval),
            chi2=chi2,
            ndof=ndof,
            x_range=self._x_range,
            at_limit=at_limit,
            minuit=m,
        )

    def plot(
        self,
        result: UnbinnedFitResult,
        hist: Any,
        npts: int = 1000,
        logy: bool = True,
        subdirectory: str = "fits",
        label: Optional[str] = None,
    ) -> None:
        """Overlay the fit model on a caller-supplied ROOT ``TH1``.

        The histogram is the source of truth for binning — the caller
        prepares it (e.g. by reading ``calibrated_zoomedHist`` from the
        ROOT file) with whatever bin width the analysis uses. The fit
        curve is scaled to match that bin width so that
        ``rate_density(x) * bin_width`` is plotted.

        Parameters
        ----------
        result : UnbinnedFitResult
        hist : ROOT.TH1
            Pre-binned data histogram to draw on. Bin width is read from
            its first bin; uniform binning is assumed. Bins outside the
            fit window are shown but ignored for residuals (handled by
            :cpp:func:`PlottingUtils::PlotFitWithResiduals`).
        npts : int, optional
            Points in the fine grid used for the overlay curve.
        logy : bool, optional
        subdirectory : str, optional
        label : str, optional
        """
        from analysis_utilities import load_cpp_library
        ROOT = load_cpp_library()

        if self._events is None or self._model is None or self._x_range is None:
            raise RuntimeError("set_events() has not been called")

        x_lo, x_hi = self._x_range
        hist_bw = float(hist.GetXaxis().GetBinWidth(1))

        fine_edges = np.linspace(x_lo, x_hi, npts + 1)
        fine_bw = (x_hi - x_lo) / npts
        fine_centers = 0.5 * (fine_edges[:-1] + fine_edges[1:])

        param_names = _model_params(self._model)
        param_values = [result.values[p] for p in param_names]
        total_y = _evaluate_binned(self._model, fine_edges, param_values,
                                    hist_bw)
        total_graph = _make_tgraph(ROOT, fine_centers, total_y, ROOT.kAzure)

        overlay_colors = [ROOT.kBlack, ROOT.kGray, ROOT.kRed, ROOT.kOrange]
        component_graphs: List[Any] = []
        if self._components:
            bkg_fine: Optional[np.ndarray] = None
            if (self._background_component is not None
                    and self._background_component in self._components):
                bkg_fn = self._components[self._background_component]
                bkg_vals = [result.values[p] for p in _model_params(bkg_fn)]
                bkg_fine = _evaluate_binned(bkg_fn, fine_edges, bkg_vals,
                                             hist_bw)

            color_cursor = 0
            for comp_name, comp_fn in self._components.items():
                comp_vals = [result.values[p] for p in _model_params(comp_fn)]
                comp_fine = _evaluate_binned(comp_fn, fine_edges, comp_vals,
                                              hist_bw)
                if comp_name == self._background_component:
                    overlay = comp_fine
                    color = ROOT.kGreen
                elif bkg_fine is not None:
                    overlay = comp_fine + bkg_fine
                    color = overlay_colors[color_cursor % len(overlay_colors)]
                    color_cursor += 1
                else:
                    overlay = comp_fine
                    color = overlay_colors[color_cursor % len(overlay_colors)]
                    color_cursor += 1
                component_graphs.append(
                    _make_tgraph(ROOT, fine_centers, overlay, color))

        comp_vec = ROOT.std.vector("TGraph*")()
        for g in component_graphs:
            comp_vec.push_back(g)

        if label is None:
            reduced = result.chi2 / max(result.ndof, 1)
            label = f"#chi^{{2}}/ndf = {reduced:.3f}"

        out_name = self.name
        ROOT.PlottingUtils.PlotFitWithResiduals(
            hist,
            total_graph,
            comp_vec,
            float(x_lo),
            float(x_hi),
            out_name,
            subdirectory,
            label,
            logy,
        )


def _model_params(model: Callable) -> List[str]:
    sig = inspect.signature(model)
    all_params = list(sig.parameters.values())
    return [p.name for p in all_params[1:]]


def _evaluate_binned(
    model: Callable,
    fine_edges: np.ndarray,
    param_values: Sequence[float],
    hist_bw: float,
) -> np.ndarray:
    """Convert ``(total, rate_density)`` model output to counts-per-hist-bin.

    Expected count at a point is ``rate_density(x) * hist_bw`` — the
    model's rate density (events-per-unit-x) times the bin width of the
    histogram being overlaid.
    """
    _, rate_density = model(fine_edges, *param_values)
    rate_at_center = 0.5 * (rate_density[:-1] + rate_density[1:])
    return rate_at_center * hist_bw


def _make_tgraph(
    ROOT: Any,
    x: Sequence[float],
    y: Sequence[float],
    color: int,
) -> Any:
    n = len(x)
    x_arr = _array.array("d", [float(v) for v in x])
    y_arr = _array.array("d", [float(v) for v in y])
    graph = ROOT.TGraph(n, x_arr, y_arr)
    graph.SetLineColor(color)
    graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    return graph
