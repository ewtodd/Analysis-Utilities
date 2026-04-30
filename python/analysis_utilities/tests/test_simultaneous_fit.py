"""Tests for analysis_utilities.simultaneous_fit."""

import numpy as np
import pytest

pytest.importorskip("iminuit")

from analysis_utilities.fit_models import (compose, exponential_bkg,
                                             gaussian_peak)
from analysis_utilities.simultaneous_fit import SimultaneousFit


def _true_sr(edges: np.ndarray, mu: float, sigma: float, nsig: float,
             tau: float, nbkg: float, ref_lo: float,
             ref_hi: float) -> np.ndarray:
    from scipy.special import erf
    z = (edges - mu) / (np.sqrt(2.0) * sigma)
    peak = nsig * np.diff(0.5 * (1.0 + erf(z)))
    shifted_hi = ref_hi - ref_lo
    ref_integral = 1.0 - np.exp(-shifted_hi / tau)
    bin_integrals = (np.exp(-(edges[:-1] - ref_lo) / tau)
                     - np.exp(-(edges[1:] - ref_lo) / tau))
    return peak + nbkg * bin_integrals / ref_integral


def _true_bkg(edges: np.ndarray, tau: float, nbkg: float, ref_lo: float,
              ref_hi: float) -> np.ndarray:
    shifted_hi = ref_hi - ref_lo
    ref_integral = 1.0 - np.exp(-shifted_hi / tau)
    bin_integrals = (np.exp(-(edges[:-1] - ref_lo) / tau)
                     - np.exp(-(edges[1:] - ref_lo) / tau))
    return nbkg * bin_integrals / ref_integral


def test_sr_sideband_recovers_truth_within_three_sigma() -> None:
    rng = np.random.default_rng(0)

    ref_lo, ref_hi = 60.0, 80.0
    edges_sr = np.linspace(60.0, 80.0, 201)
    edges_sb = np.linspace(85.0, 100.0, 151)

    mu_true, sigma_true, nsig_true = 68.75, 0.4, 500.0
    tau_true = 40.0
    nbkg_sr_true, nbkg_sb_true = 10000.0, 3000.0

    counts_sr = rng.poisson(
        _true_sr(edges_sr, mu_true, sigma_true, nsig_true, tau_true,
                 nbkg_sr_true, ref_lo, ref_hi)).astype(np.float64)
    counts_sb = rng.poisson(
        _true_bkg(edges_sb, tau_true, nbkg_sb_true, ref_lo,
                  ref_hi)).astype(np.float64)

    peak = gaussian_peak(mu="mu", sigma="sigma", nsig="nsig")
    bkg_sr = exponential_bkg((ref_lo, ref_hi), tau="tau", nbkg="nbkg_sr")
    bkg_sb = exponential_bkg((ref_lo, ref_hi), tau="tau", nbkg="nbkg_sb")

    sr_model = compose(signal=peak, background=bkg_sr)

    fit = SimultaneousFit(name="test_srsb")
    fit.add_channel("sr",
                    counts_sr,
                    edges_sr,
                    sr_model,
                    background_component="background")
    fit.add_channel("sb", counts_sb, edges_sb, bkg_sb)

    fit.set_parameter("mu", value=68.5, limits=(60.0, 80.0))
    fit.set_parameter("sigma", value=0.5, limits=(0.05, 5.0))
    fit.set_parameter("nsig", value=400.0, limits=(0.0, None))
    fit.set_parameter("tau", value=30.0, limits=(1.0, 500.0))
    fit.set_parameter("nbkg_sr", value=9000.0, limits=(0.0, None))
    fit.set_parameter("nbkg_sb", value=3000.0, limits=(0.0, None))

    result = fit.fit()

    assert result.valid
    truth = {
        "mu": mu_true,
        "sigma": sigma_true,
        "nsig": nsig_true,
        "tau": tau_true,
        "nbkg_sr": nbkg_sr_true,
        "nbkg_sb": nbkg_sb_true,
    }
    for name, true_val in truth.items():
        pull = (result.values[name] - true_val) / result.errors[name]
        assert abs(pull) < 3.0, f"{name}: pull {pull:+.2f}"


def test_shared_parameter_counted_once() -> None:
    rng = np.random.default_rng(1)

    ref = (0.0, 10.0)
    edges_a = np.linspace(0.0, 10.0, 51)
    edges_b = np.linspace(0.0, 10.0, 51)

    shape_a = exponential_bkg(ref, tau="tau_shared", nbkg="n_a")
    shape_b = exponential_bkg(ref, tau="tau_shared", nbkg="n_b")

    counts_a = rng.poisson(shape_a(edges_a, 3.0, 1000.0)).astype(np.float64)
    counts_b = rng.poisson(shape_b(edges_b, 3.0, 500.0)).astype(np.float64)

    fit = SimultaneousFit(name="test_share")
    fit.add_channel("a", counts_a, edges_a, shape_a)
    fit.add_channel("b", counts_b, edges_b, shape_b)
    fit.set_parameter("tau_shared", value=2.5, limits=(0.1, 50.0))
    fit.set_parameter("n_a", value=900.0, limits=(0.0, None))
    fit.set_parameter("n_b", value=450.0, limits=(0.0, None))

    result = fit.fit()
    assert result.valid
    assert set(result.values.keys()) == {"tau_shared", "n_a", "n_b"}


def test_normal_constraint_pulls_parameter() -> None:
    rng = np.random.default_rng(2)
    ref = (0.0, 10.0)
    edges = np.linspace(0.0, 10.0, 51)
    shape = exponential_bkg(ref, tau="tau", nbkg="nbkg")
    counts = rng.poisson(shape(edges, 3.0, 1000.0)).astype(np.float64)

    fit = SimultaneousFit(name="test_cons")
    fit.add_channel("a", counts, edges, shape)
    fit.set_parameter("tau", value=2.5, limits=(0.1, 50.0))
    fit.set_parameter("nbkg", value=900.0, limits=(0.0, None))
    fit.add_constraint("tau", mean=5.0, sigma=0.1)

    result = fit.fit()
    assert result.valid
    assert abs(result.values["tau"] - 5.0) < 0.5
