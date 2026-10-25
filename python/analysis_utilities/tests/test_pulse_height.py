"""PulseHeightProcessor's Deconvolve and Trapezoid against their definitions.

Each filter is checked against a closed form where one exists, and against an
independent direct computation otherwise, so the formulas are verified
explicitly rather than through the figures they produce.
"""

from __future__ import annotations

import math
import random

import pytest

ROOT = pytest.importorskip("ROOT")

try:
    from analysis_utilities import load_cpp_library

    load_cpp_library()
except (RuntimeError, OSError) as exc:
    pytest.skip(f"C++ library not loadable: {exc}", allow_module_level=True)


def _to_tarray(values: list[float]):
    array = ROOT.TArrayF(len(values))
    for index, value in enumerate(values):
        array.AddAt(value, index)
    return array


def _to_list(array) -> list[float]:
    return [array.At(index) for index in range(array.GetSize())]


def _trapezoid(values: list[float], rise: int, flat: int) -> list[float]:
    return _to_list(ROOT.PulseHeightProcessor.Trapezoid(
        _to_tarray(values), rise, flat))


def _deconvolve(values: list[float], tau: float) -> list[float]:
    return _to_list(ROOT.PulseHeightProcessor.Deconvolve(
        _to_tarray(values), tau))


def test_trapezoid_of_a_step_matches_the_closed_form() -> None:
    """A step of height A at sample s gives, with samples before the start at 0,
    T[n] = A (clamp(n - s + 1) - clamp(n - k - m - s + 1)) / k, clamp to [0, k].
    """
    n_samples, start, height, rise, flat = 700, 60, 1000.0, 100, 40
    step = [height if n >= start else 0.0 for n in range(n_samples)]
    trapezoid = _trapezoid(step, rise, flat)

    def clamp(x: int) -> int:
        return min(max(x, 0), rise)

    for n in range(n_samples):
        expected = height * (clamp(n - start + 1)
                             - clamp(n - rise - flat - start + 1)) / rise
        assert trapezoid[n] == pytest.approx(expected, abs=1e-3), f"sample {n}"


def test_trapezoid_flat_top_spans_s_plus_k_minus_1_to_s_plus_k_plus_m_minus_1(
) -> None:
    n_samples, start, height, rise, flat = 700, 60, 1000.0, 100, 40
    step = [height if n >= start else 0.0 for n in range(n_samples)]
    trapezoid = _trapezoid(step, rise, flat)
    first, last = start + rise - 1, start + rise + flat - 1
    for n in range(first, last + 1):
        assert trapezoid[n] == pytest.approx(height, abs=1e-3), f"sample {n}"
    assert trapezoid[first - 1] < height - 1.0
    assert trapezoid[last + 1] < height - 1.0


def test_trapezoid_matches_a_direct_sum() -> None:
    """Random input against the definition summed sample by sample: the mean
    over the last k samples minus the mean over the k samples ending k + m
    earlier, with samples before the start counted as zero."""
    rng = random.Random(20261024)
    n_samples, rise, flat = 500, 37, 13
    values = [rng.gauss(0.0, 30.0) + (400.0 if n > 200 else 0.0)
              for n in range(n_samples)]
    trapezoid = _trapezoid(values, rise, flat)
    for n in range(n_samples):
        late = sum(values[j] for j in range(max(0, n - rise + 1), n + 1))
        early_last = n - rise - flat
        early = sum(values[j]
                    for j in range(max(0, early_last - rise + 1),
                                   early_last + 1))
        expected = (late - early) / rise
        assert trapezoid[n] == pytest.approx(expected, rel=1e-5, abs=1e-3), (
            f"sample {n}")


def test_deconvolve_inverts_a_discrete_exponential() -> None:
    """x[n] = d x[n-1] + i[n], d = exp(-1/tau), deconvolves to sum of i."""
    rng = random.Random(7)
    n_samples, tau = 800, 324.0
    decay = math.exp(-1.0 / tau)
    impulses = [rng.uniform(50.0, 500.0) if rng.random() < 0.01 else 0.0
                for n in range(n_samples)]
    signal, previous = [], 0.0
    for impulse in impulses:
        previous = decay * previous + impulse
        signal.append(previous)
    step = _deconvolve(signal, tau)
    running = 0.0
    for n in range(n_samples):
        running += impulses[n]
        assert step[n] == pytest.approx(running, rel=1e-4, abs=1e-2), (
            f"sample {n}")


def test_deconvolve_without_a_decay_constant_returns_the_input() -> None:
    values = [float(n % 17) for n in range(100)]
    assert _deconvolve(values, 0.0) == pytest.approx(values)


def test_trapezoid_of_a_deconvolved_pulse_reads_its_amplitude() -> None:
    n_samples, start, amplitude, tau, rise, flat = 1500, 300, 5000.0, 324.0, \
        400, 100
    decay = math.exp(-1.0 / tau)
    pulse = [amplitude * decay ** (n - start) if n >= start else 0.0
             for n in range(n_samples)]
    trapezoid = _trapezoid(_deconvolve(pulse, tau), rise, flat)
    for n in range(start + rise - 1, start + rise + flat):
        assert trapezoid[n] == pytest.approx(amplitude, rel=1e-5), (
            f"sample {n}")
