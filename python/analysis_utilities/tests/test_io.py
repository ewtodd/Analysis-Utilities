"""Tests for analysis_utilities.io read/write helpers."""

from __future__ import annotations

import warnings

import pytest

ROOT = pytest.importorskip("ROOT")

try:
    from analysis_utilities import (get_root_files_base_dir, load_cpp_library,
                                      open_for_reading, open_for_writing,
                                      set_root_files_base_dir,
                                      set_root_preferences)

    load_cpp_library()
except (RuntimeError, OSError) as exc:
    pytest.skip(f"C++ library not loadable: {exc}", allow_module_level=True)


def test_open_for_writing_creates_nested_parents(tmp_path) -> None:
    set_root_files_base_dir(tmp_path)
    fout = open_for_writing("deep/nested/dir/file.root")
    try:
        assert fout.IsOpen()
    finally:
        fout.Close()
    expected = tmp_path / "deep" / "nested" / "dir" / "file.root"
    assert expected.exists()


def test_round_trip_read_write(tmp_path) -> None:
    set_root_files_base_dir(tmp_path)

    fout = open_for_writing("round/trip.root")
    h = ROOT.TH1F("h", "h", 10, 0.0, 1.0)
    h.SetBinContent(5, 42.0)
    h.Write()
    fout.Close()

    fin = open_for_reading("round/trip.root")
    try:
        assert fin.IsOpen()
        h_read = fin.Get("h")
        assert h_read is not None
        assert h_read.GetBinContent(5) == 42.0
    finally:
        fin.Close()


def test_absolute_path_passthrough(tmp_path) -> None:
    other_base = tmp_path / "ignored_base"
    set_root_files_base_dir(other_base)

    abs_target = tmp_path / "abs_target.root"
    fout = open_for_writing(str(abs_target))
    try:
        assert fout.IsOpen()
    finally:
        fout.Close()
    assert abs_target.exists()
    # Absolute path should bypass the base dir entirely; base dir was
    # never touched, so it shouldn't have been created.
    assert not other_base.exists()


def test_set_root_preferences_propagates_root_files_base(tmp_path) -> None:
    plots_dir = tmp_path / "plots"
    root_files_dir = tmp_path / "root_files"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        set_root_preferences(plots_dir=plots_dir, root_files_dir=root_files_dir)
    assert get_root_files_base_dir() == str(root_files_dir.resolve())
    assert (
        str(ROOT.IO.GetRootFilesBaseDir()) == str(root_files_dir.resolve())
    )


def test_set_root_files_base_dir_strips_trailing_slash(tmp_path) -> None:
    target = tmp_path / "no_slash"
    set_root_files_base_dir(str(target) + "/")
    assert get_root_files_base_dir() == str(target.resolve())
    assert str(ROOT.IO.GetRootFilesBaseDir()) == str(target.resolve())


def _write_leaf_array_tree(path, n_events) -> None:
    """Write a tree with two fixed-size leaf-list arrays and a scalar."""
    import array as _array

    fout = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("events", "leaf-array test tree")
    left = _array.array("H", [0] * 4)
    right = _array.array("f", [0.0] * 4)
    cathode = _array.array("h", [0])
    tree.Branch("Left", left, "Left[4]/s")
    tree.Branch("Right", right, "Right[4]/F")
    tree.Branch("Cathode", cathode, "Cathode/S")
    for i in range(n_events):
        for s in range(4):
            left[s] = 10 * i + s
            right[s] = 0.5 * (10 * i + s)
        cathode[0] = -i
        tree.Fill()
    tree.Write()
    fout.Close()


def test_load_leaf_array_data_round_trip(tmp_path) -> None:
    import numpy as np

    from analysis_utilities.io import load_leaf_array_data

    root_path = tmp_path / "leafarr.root"
    _write_leaf_array_tree(root_path, n_events=5)

    arrays = load_leaf_array_data(str(root_path), "events",
                                  ["Left", "Right"],
                                  cache_dir=str(tmp_path / "cache"))
    assert set(arrays) == {"Left", "Right"}
    assert arrays["Left"].shape == (5, 4)
    assert arrays["Left"].dtype == np.uint16
    assert arrays["Right"].dtype == np.float32
    assert arrays["Left"][3, 2] == 32
    assert arrays["Right"][3, 2] == pytest.approx(16.0)

    # Second call must hit the .npz cache and return identical contents.
    cached = load_leaf_array_data(str(root_path), "events",
                                  ["Left", "Right"],
                                  cache_dir=str(tmp_path / "cache"))
    assert np.array_equal(cached["Left"], arrays["Left"])
    assert np.array_equal(cached["Right"], arrays["Right"])


def test_load_leaf_array_data_rejects_scalar_branch(tmp_path) -> None:
    from analysis_utilities.io import load_leaf_array_data

    root_path = tmp_path / "leafarr_scalar.root"
    _write_leaf_array_tree(root_path, n_events=2)

    with pytest.raises(ValueError, match="not a fixed-size array"):
        load_leaf_array_data(str(root_path), "events", ["Cathode"],
                             cache_dir=None)


def test_load_leaf_array_data_max_events(tmp_path) -> None:
    from analysis_utilities.io import load_leaf_array_data

    root_path = tmp_path / "leafarr_cap.root"
    _write_leaf_array_tree(root_path, n_events=6)

    arrays = load_leaf_array_data(str(root_path), "events", ["Left"],
                                  max_events=4, cache_dir=None)
    assert arrays["Left"].shape == (4, 4)
