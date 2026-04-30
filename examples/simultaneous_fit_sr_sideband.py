"""Simultaneous fit example: Ge peak on Pb-K-alpha shoulder.

Two channels share an exponential background shape:

  * ``signal_region``: gaussian peak at ~68.75 keV on exponential background.
  * ``sideband``: exponential-only window at higher energy that constrains
    the shared decay constant ``tau``.

The sideband has no peak parameters, so the only coupling between channels
is the shared ``tau``. This is the canonical case where a simultaneous fit
buys you something a single-channel fit does not: the sideband data
constrains ``tau``, which in turn reduces the ``tau`` <-> ``nbkg_sr``
correlation in the signal region and tightens the ``nsig`` uncertainty.

Run inside a ``nix develop`` shell (or any environment with the
``analysis-utilities`` package + ROOT available).
"""

from __future__ import annotations

import numpy as np

from analysis_utilities import SimultaneousFit, load_cpp_library
from analysis_utilities.fit_models import (compose, exponential_bkg,
                                             gaussian_peak)


def generate_toy_data(seed: int = 123):
    rng = np.random.default_rng(seed)

    ref_lo, ref_hi = 60.0, 80.0

    edges_sr = np.linspace(65.0, 75.0, 201)
    edges_sb = np.linspace(85.0, 100.0, 151)

    mu_true, sigma_true, nsig_true = 68.75, 0.35, 600.0
    tau_true = 50.0
    nbkg_sr_true, nbkg_sb_true = 12000.0, 4000.0

    from scipy.special import erf

    def peak_counts(edges: np.ndarray) -> np.ndarray:
        z = (edges - mu_true) / (np.sqrt(2.0) * sigma_true)
        return nsig_true * np.diff(0.5 * (1.0 + erf(z)))

    def bkg_counts(edges: np.ndarray, nbkg: float) -> np.ndarray:
        shifted_hi = ref_hi - ref_lo
        ref_integral = 1.0 - np.exp(-shifted_hi / tau_true)
        bin_integrals = (np.exp(-(edges[:-1] - ref_lo) / tau_true)
                         - np.exp(-(edges[1:] - ref_lo) / tau_true))
        return nbkg * bin_integrals / ref_integral

    expected_sr = peak_counts(edges_sr) + bkg_counts(edges_sr, nbkg_sr_true)
    expected_sb = bkg_counts(edges_sb, nbkg_sb_true)

    counts_sr = rng.poisson(expected_sr).astype(np.float64)
    counts_sb = rng.poisson(expected_sb).astype(np.float64)

    return (edges_sr, counts_sr), (edges_sb, counts_sb), {
        "mu": mu_true,
        "sigma": sigma_true,
        "nsig": nsig_true,
        "tau": tau_true,
        "nbkg_sr": nbkg_sr_true,
        "nbkg_sb": nbkg_sb_true,
        "reference_range": (ref_lo, ref_hi),
    }


def main() -> None:
    ROOT = load_cpp_library()
    ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)

    (edges_sr, counts_sr), (edges_sb, counts_sb), truth = generate_toy_data()
    ref_range = truth["reference_range"]

    peak = gaussian_peak(mu="mu", sigma="sigma", nsig="nsig")
    bkg_sr = exponential_bkg(ref_range, tau="tau", nbkg="nbkg_sr")
    bkg_sb = exponential_bkg(ref_range, tau="tau", nbkg="nbkg_sb")
    sr_model = compose(signal=peak, background=bkg_sr)

    fit = SimultaneousFit(name="ge73m_demo")
    fit.add_channel("signal_region",
                    counts_sr,
                    edges_sr,
                    sr_model,
                    background_component="background")
    fit.add_channel("sideband", counts_sb, edges_sb, bkg_sb)

    fit.set_parameter("mu", value=68.5, limits=(65.0, 75.0))
    fit.set_parameter("sigma", value=0.5, limits=(0.05, 5.0))
    fit.set_parameter("nsig", value=400.0, limits=(0.0, None))
    fit.set_parameter("tau", value=30.0, limits=(1.0, 500.0))
    fit.set_parameter("nbkg_sr", value=10000.0, limits=(0.0, None))
    fit.set_parameter("nbkg_sb", value=3500.0, limits=(0.0, None))

    result = fit.fit(minos=["nsig", "mu"])

    print(f"valid={result.valid}  "
          f"nll={result.nll:.3f}  "
          f"chi2/ndof = {result.chi2:.1f} / {result.ndof}")
    for name in ("mu", "sigma", "nsig", "tau", "nbkg_sr", "nbkg_sb"):
        v = result.values[name]
        e = result.errors[name]
        true_val = truth[name]
        pull = (v - true_val) / e
        extra = ""
        if name in result.minos:
            lo, hi = result.minos[name]
            extra = f"  [minos {lo:+.3f} / {hi:+.3f}]"
        print(f"  {name:10s} = {v:10.4f} +/- {e:.4f}   "
              f"(truth {true_val:.4f}, pull {pull:+.2f}){extra}")

    fit.plot(
        result,
        labels={
            "signal_region": "signal region",
            "sideband": "sideband",
        },
    )


if __name__ == "__main__":
    main()
