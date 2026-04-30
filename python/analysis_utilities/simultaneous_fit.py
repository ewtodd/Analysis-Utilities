"""Simultaneous binned fit across channels with name-based parameter sharing.

Build combined extended binned likelihoods from per-channel parametric
models, add optional Gaussian constraints, run Minuit, and produce
publication-styled plots via :cpp:class:`PlottingUtils::PlotFitWithResiduals`.

Complements :cpp:class:`FittingUtils` (single-channel peak fits with
component pruning) for problems with shared parameters across channels —
e.g. signal region + sideband with a shared background shape.
"""

from __future__ import annotations

import array as _array
import inspect
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass
class _ChannelSpec:
    name: str
    counts: np.ndarray
    edges: np.ndarray
    model: Callable
    components: Optional[Dict[str, Callable]]
    background_component: Optional[str]
    th1: Any


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
class FitResult:
    """Structured fit output.

    Attributes
    ----------
    valid : bool
        Minuit reported a valid minimum and a non-null covariance.
    values : dict[str, float]
        Central parameter values at the minimum.
    errors : dict[str, float]
        Symmetric (Hesse) parameter errors.
    minos : dict[str, tuple[float, float]]
        Asymmetric (Minos) errors as ``(lower, upper)``, only for parameters
        for which Minos was requested. Empty if ``minos`` was not passed to
        :meth:`SimultaneousFit.fit`.
    corr : pandas.DataFrame
        Correlation matrix of free parameters, indexed by parameter name.
    nll : float
        Value of the cost function (2 * negative log likelihood) at the
        minimum, summed across channels and constraints.
    chi2 : float
        Pearson chi-squared summed across channels:
        :math:`\\sum_i (N_i - \\hat{E}_i)^2 / \\hat{E}_i` over all bins.
    ndof : int
        Total bin count minus number of free parameters.
    chi2_per_channel : dict[str, tuple[float, int]]
        Per-channel ``(chi2, n_bins)``. Degrees of freedom are not
        attributed per channel because parameters are shared.
    residuals : dict[str, tuple[ndarray, ndarray, ndarray]]
        Per-channel ``(bin_centers, pulls, raw_residuals)`` where
        ``pulls = (N - E) / sqrt(E)`` and ``raw_residuals = N - E``.
    at_limit : list[str]
        Free parameters whose fitted values are within ``1e-4`` (relative)
        of a limit. Each entry is ``"<name>@lower=<val>"`` or
        ``"<name>@upper=<val>"``. Empty if no parameters are at their
        limits. A common cause of Minuit invalidity.
    minuit : iminuit.Minuit
        The underlying Minuit instance (escape hatch for ``mnprofile``,
        ``draw_mnmatrix``, etc.).
    """

    valid: bool
    values: Dict[str, float]
    errors: Dict[str, float]
    minos: Dict[str, Tuple[float, float]]
    corr: pd.DataFrame
    nll: float
    chi2: float
    ndof: int
    chi2_per_channel: Dict[str, Tuple[float, int]]
    residuals: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]
    at_limit: List[str]
    minuit: Any


class SimultaneousFit:
    """Combined extended binned likelihood fit over multiple channels.

    Parameters sharing between channels is by name: any parameter appearing
    in multiple channel model signatures is a single fit parameter. This is
    native :mod:`iminuit.cost` behavior — the combined cost is built with
    the ``+`` operator.

    Parameters
    ----------
    name : str, optional
        Output-name prefix for plots. Per-channel files are written to
        ``plots/<subdirectory>/<name>_<channel>.{png|pdf}``.

    See Also
    --------
    :class:`FittingUtils` (C++) : Single-channel peak fits with component
        pruning (step, exponential tails) and interactive editor.
    """

    def __init__(self, name: str = "simfit") -> None:
        self.name: str = name
        self._channels: List[_ChannelSpec] = []
        self._param_specs: Dict[str, _ParameterSpec] = {}
        self._constraints: List[_ConstraintSpec] = []

    def add_channel(
        self,
        name: str,
        counts: np.ndarray,
        edges: np.ndarray,
        model: Callable,
        components: Optional[Dict[str, Callable]] = None,
        background_component: Optional[str] = None,
    ) -> None:
        """Add a channel from raw counts and bin edges.

        Parameters
        ----------
        name : str
            Unique channel identifier.
        counts : array-like
            Bin contents, length ``N``.
        edges : array-like
            Bin edges, length ``N+1``. Uniform spacing is required for
            correct plot scaling.
        model : callable
            ``model(edges, *params) -> ndarray`` returning per-bin expected
            counts of length ``N``. Parameter names from the callable's
            signature are used for parameter sharing. Builders from
            :mod:`fit_models` produce such callables.
        components : dict[str, callable], optional
            Named sub-components used as overlays in plots. If omitted and
            ``model`` has a ``.components`` attribute (as produced by
            :func:`fit_models.compose`), that is used.
        background_component : str, optional
            Key within ``components`` identifying the background. Non-
            background components are drawn summed with the background,
            matching the :cpp:class:`FittingUtils` convention.
        """
        counts_arr = np.asarray(counts, dtype=np.float64)
        edges_arr = np.asarray(edges, dtype=np.float64)
        if len(edges_arr) != len(counts_arr) + 1:
            raise ValueError(
                f"edges length must be counts length + 1 "
                f"(got {len(edges_arr)} vs {len(counts_arr)})")
        if components is None:
            components = getattr(model, "components", None)
        self._channels.append(
            _ChannelSpec(
                name=name,
                counts=counts_arr,
                edges=edges_arr,
                model=model,
                components=components,
                background_component=background_component,
                th1=None,
            ))

    def add_channel_from_th1(
        self,
        name: str,
        hist: Any,
        model: Callable,
        x_range: Optional[Tuple[float, float]] = None,
        components: Optional[Dict[str, Callable]] = None,
        background_component: Optional[str] = None,
    ) -> None:
        """Add a channel from a ROOT ``TH1``.

        Parameters
        ----------
        name : str
        hist : ROOT.TH1
            1-D ROOT histogram. A reference is kept for plotting.
        model : callable
        x_range : tuple(float, float), optional
            Restrict to bins fully contained in ``(lo, hi)``.
        components, background_component
            See :meth:`add_channel`.
        """
        nbins = int(hist.GetNbinsX())
        axis = hist.GetXaxis()
        edges = np.empty(nbins + 1, dtype=np.float64)
        for i in range(nbins):
            edges[i] = float(axis.GetBinLowEdge(i + 1))
        edges[nbins] = float(axis.GetBinUpEdge(nbins))
        counts = np.empty(nbins, dtype=np.float64)
        for i in range(nbins):
            counts[i] = float(hist.GetBinContent(i + 1))

        if x_range is not None:
            lo, hi = float(x_range[0]), float(x_range[1])
            eps = 1e-9 * (hi - lo)
            mask = (edges[:-1] >= lo - eps) & (edges[1:] <= hi + eps)
            if not np.any(mask):
                raise ValueError(
                    f"x_range ({lo}, {hi}) selects no bins from {name}")
            first = int(np.argmax(mask))
            last = len(mask) - 1 - int(np.argmax(mask[::-1]))
            counts = counts[first:last + 1].copy()
            edges = edges[first:last + 2].copy()

        if components is None:
            components = getattr(model, "components", None)

        self._channels.append(
            _ChannelSpec(
                name=name,
                counts=counts,
                edges=edges,
                model=model,
                components=components,
                background_component=background_component,
                th1=hist,
            ))

    def set_parameter(
        self,
        name: str,
        value: Optional[float] = None,
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        fixed: bool = False,
        error: Optional[float] = None,
    ) -> None:
        """Set initial value, limits, fixed status, and step for a parameter.

        Parameters
        ----------
        name : str
            Parameter name as it appears in one or more channel model
            signatures.
        value : float, optional
            Initial value. If omitted, defaults to 1.0 at fit time (a
            warning is emitted).
        limits : tuple, optional
            ``(lower, upper)``; either may be ``None`` for unbounded.
        fixed : bool, optional
            Fix the parameter at ``value`` during the fit.
        error : float, optional
            Initial step size (MINUIT "error"). Defaults to iminuit's
            automatic choice.
        """
        self._param_specs[name] = _ParameterSpec(
            value=value,
            limits=limits,
            fixed=fixed,
            error=error,
        )

    def add_constraint(self, name: str, mean: float, sigma: float) -> None:
        """Add a Gaussian constraint term to the cost function.

        Implemented via :class:`iminuit.cost.NormalConstraint`, added to the
        combined cost with ``+``.
        """
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
    ) -> FitResult:
        """Run ``migrad`` + ``hesse`` (+ optional ``minos``).

        Parameters
        ----------
        minos : bool or iterable of str, optional
            If ``True``, run Minos on every free parameter. If an iterable
            of names, run Minos only on those. If ``None`` (default), skip.
        strategy : int, optional
            Minuit strategy, 0 (fast), 1 (default), or 2 (careful).
        tol : float, optional
            Minuit tolerance. Uses iminuit's default if omitted.
        print_level : int, optional
            Minuit verbosity (0 silent).
        migrad_ncall : int, optional
            Max function calls per migrad pass. ``None`` (default) uses
            iminuit's automatic heuristic.
        migrad_iterate : int, optional
            Restart migrad up to this many times while :attr:`Minuit.valid`
            is ``False`` but further progress is being made. Passed to
            :meth:`iminuit.Minuit.migrad`.

        Returns
        -------
        FitResult
        """
        from iminuit import Minuit
        from iminuit.cost import ExtendedBinnedNLL, NormalConstraint

        if not self._channels:
            raise RuntimeError("fit() called with no channels")

        combined_cost: Any = None
        for ch in self._channels:
            scaled_cdf = _make_scaled_cdf(ch.model)
            channel_cost = ExtendedBinnedNLL(ch.counts, ch.edges, scaled_cdf)
            combined_cost = (channel_cost if combined_cost is None
                             else combined_cost + channel_cost)

        for con in self._constraints:
            nc = NormalConstraint(con.name, con.mean, con.sigma)
            combined_cost = nc if combined_cost is None else combined_cost + nc

        all_params: List[str] = list(combined_cost._parameters.keys())

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

        m = Minuit(combined_cost, **initial_values)
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

        chi2_per_channel: Dict[str, Tuple[float, int]] = {}
        residuals: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        chi2_total: float = 0.0
        total_bins: int = 0

        for ch in self._channels:
            param_names = _model_params(ch.model)
            param_values = [values[p] for p in param_names]
            expected = np.asarray(ch.model(ch.edges, *param_values),
                                   dtype=np.float64)
            centers = 0.5 * (ch.edges[:-1] + ch.edges[1:])
            with np.errstate(divide="ignore", invalid="ignore"):
                pulls = np.where(
                    expected > 0,
                    (ch.counts - expected) / np.sqrt(expected),
                    0.0,
                )
            ch_chi2 = float(np.sum(pulls * pulls))
            chi2_per_channel[ch.name] = (ch_chi2, int(len(ch.counts)))
            chi2_total += ch_chi2
            total_bins += len(ch.counts)
            residuals[ch.name] = (centers, pulls,
                                  (ch.counts - expected))

        ndof_total: int = total_bins - len(free_params)

        return FitResult(
            valid=bool(m.valid and m.covariance is not None),
            values=values,
            errors=errors,
            minos=minos_results,
            corr=corr,
            nll=float(m.fval),
            chi2=chi2_total,
            ndof=ndof_total,
            chi2_per_channel=chi2_per_channel,
            residuals=residuals,
            at_limit=at_limit,
            minuit=m,
        )

    def plot(
        self,
        result: FitResult,
        npts: int = 1000,
        logy: bool = True,
        subdirectory: str = "fits",
        labels: Optional[Dict[str, str]] = None,
    ) -> None:
        """Produce a per-channel fit plot via PlottingUtils.

        One two-pad canvas per channel (main + residual) plus a pull
        histogram, written to ``plots/<subdirectory>/``.

        Parameters
        ----------
        result : FitResult
        npts : int, optional
            Points in the fine grid used for overlay curves.
        logy : bool, optional
            Log-y on the main pad.
        subdirectory : str, optional
            Passed through to :cpp:func:`PlottingUtils::SaveFigure`.
        labels : dict[str, str], optional
            Per-channel annotation label drawn on the main pad. If
            omitted, every panel is labeled with the global reduced
            chi-squared ``"#chi^{2}/ndf = %.3f"``, matching the
            :cpp:class:`FittingUtils` convention.

        Notes
        -----
        Requires :func:`analysis_utilities.load_cpp_library`. Assumes
        uniform bin widths per channel; variable-width binning will produce
        a correct likelihood but a miscalibrated plot overlay.
        """
        from analysis_utilities import load_cpp_library
        ROOT = load_cpp_library()

        reduced_chi2 = result.chi2 / max(result.ndof, 1)
        default_label = f"#chi^{{2}}/ndf = {reduced_chi2:.3f}"
        if labels is None:
            labels = {ch.name: default_label for ch in self._channels}

        overlay_colors = [ROOT.kBlack, ROOT.kGray, ROOT.kRed, ROOT.kOrange]

        for ch in self._channels:
            hist = self._build_th1(ROOT, ch)

            x_lo = float(ch.edges[0])
            x_hi = float(ch.edges[-1])
            fine_edges = np.linspace(x_lo, x_hi, npts + 1)
            fine_bw = (x_hi - x_lo) / npts
            fine_centers = 0.5 * (fine_edges[:-1] + fine_edges[1:])

            hist_bw = float(ch.edges[1] - ch.edges[0])
            scale = hist_bw / fine_bw

            total_param_names = _model_params(ch.model)
            total_param_vals = [result.values[p] for p in total_param_names]
            total_fine = np.asarray(ch.model(fine_edges, *total_param_vals),
                                     dtype=np.float64)
            total_graph = _make_tgraph(
                ROOT,
                fine_centers,
                total_fine * scale,
                color=ROOT.kAzure,
            )

            component_graphs: List[Any] = []
            if ch.components:
                bkg_fine: Optional[np.ndarray] = None
                if (ch.background_component is not None
                        and ch.background_component in ch.components):
                    bkg_fn = ch.components[ch.background_component]
                    bkg_param_names = _model_params(bkg_fn)
                    bkg_param_vals = [result.values[p]
                                       for p in bkg_param_names]
                    bkg_fine = np.asarray(
                        bkg_fn(fine_edges, *bkg_param_vals),
                        dtype=np.float64,
                    )

                color_cursor = 0
                for comp_name, comp_fn in ch.components.items():
                    comp_param_names = _model_params(comp_fn)
                    comp_param_vals = [result.values[p]
                                        for p in comp_param_names]
                    comp_fine = np.asarray(
                        comp_fn(fine_edges, *comp_param_vals),
                        dtype=np.float64,
                    )
                    if comp_name == ch.background_component:
                        overlay = comp_fine
                        color = ROOT.kGreen
                    elif bkg_fine is not None:
                        overlay = comp_fine + bkg_fine
                        color = overlay_colors[color_cursor
                                               % len(overlay_colors)]
                        color_cursor += 1
                    else:
                        overlay = comp_fine
                        color = overlay_colors[color_cursor
                                               % len(overlay_colors)]
                        color_cursor += 1
                    g = _make_tgraph(
                        ROOT,
                        fine_centers,
                        overlay * scale,
                        color=color,
                    )
                    component_graphs.append(g)

            comp_vec = ROOT.std.vector("TGraph*")()
            for g in component_graphs:
                comp_vec.push_back(g)

            out_name = f"{self.name}_{ch.name}"
            label = labels.get(ch.name, "")
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

    def _build_th1(self, ROOT: Any, ch: _ChannelSpec) -> Any:
        if ch.th1 is not None:
            return ch.th1
        edges_arr = _array.array("d", [float(v) for v in ch.edges])
        hist_name = f"_simfit_{self.name}_{ch.name}"
        n_bins = len(ch.counts)
        hist = ROOT.TH1D(hist_name, ";;Counts", n_bins, edges_arr)
        for i in range(n_bins):
            content = float(ch.counts[i])
            hist.SetBinContent(i + 1, content)
            hist.SetBinError(i + 1,
                              float(np.sqrt(max(content, 0.0))))
        return hist


def _model_params(model: Callable) -> List[str]:
    sig = inspect.signature(model)
    all_params = list(sig.parameters.values())
    return [p.name for p in all_params[1:]]


def _make_scaled_cdf(per_bin_model: Callable) -> Callable:
    signature = inspect.signature(per_bin_model)

    def wrapped(edges: np.ndarray, *params: float) -> np.ndarray:
        per_bin = per_bin_model(edges, *params)
        out = np.empty(len(edges), dtype=np.float64)
        out[0] = 0.0
        out[1:] = np.cumsum(per_bin)
        return out

    wrapped.__signature__ = signature
    return wrapped


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
