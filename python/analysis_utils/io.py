"""Generic TTree to numpy/pandas loader for waveform analysis."""

import array as _array

import numpy as np
import pandas as pd
import ROOT

# ROOT type name -> (array.array typecode, numpy dtype)
_TYPE_MAP = {
    "Float_t": ("f", np.float32),
    "float": ("f", np.float32),
    "Double_t": ("d", np.float64),
    "double": ("d", np.float64),
    "Long64_t": ("q", np.int64),
    "long long": ("q", np.int64),
    "ULong64_t": ("Q", np.uint64),
    "unsigned long long": ("Q", np.uint64),
    "Int_t": ("i", np.int32),
    "int": ("i", np.int32),
    "UInt_t": ("I", np.uint32),
    "unsigned int": ("I", np.uint32),
    "Short_t": ("h", np.int16),
    "short": ("h", np.int16),
    "UShort_t": ("H", np.uint16),
    "unsigned short": ("H", np.uint16),
    "UChar_t": ("B", np.uint8),
    "unsigned char": ("B", np.uint8),
}


def load_tree_data(
    root_files,
    tree_name="features",
    scalar_branches=None,
    array_branch=None,
    max_events=None,
):
    """Load TTree data into numpy arrays and pandas DataFrame.

    Parameters
    ----------
    root_files : str or list of str
        Path(s) to ROOT file(s).
    tree_name : str
        Name of the TTree to read.
    scalar_branches : list of str or None
        Scalar branch names to load. If None, auto-detects all
        scalar (non-array) branches.
    array_branch : str or None
        Name of a TArrayF/TArrayS branch to load as a 2-D numpy array.
    max_events : int or None
        Maximum number of events to load.

    Returns
    -------
    features_df : pandas.DataFrame
        Scalar branch data.
    waveforms : numpy.ndarray or None
        2-D array (n_events, n_samples) if array_branch is given, else None.
        Only returned if array_branch is not None.
    """
    if isinstance(root_files, str):
        root_files = [root_files]

    chain = ROOT.TChain(tree_name)
    for path in root_files:
        if chain.Add(path) == 0:
            raise FileNotFoundError(f"Could not add {path} to TChain")

    n_total = chain.GetEntries()
    if n_total == 0:
        raise ValueError(f"TChain is empty (files: {root_files})")

    # Auto-detect scalar branches if not specified
    if scalar_branches is None:
        scalar_branches = []
        branch_list = chain.GetListOfBranches()
        for i in range(branch_list.GetEntries()):
            br = branch_list.At(i)
            name = br.GetName()
            if name != array_branch:
                scalar_branches.append(name)

    if max_events is not None:
        n_to_read = min(n_total, max_events)
    else:
        n_to_read = n_total

    # Disable all branches, then enable only the ones we need
    chain.SetBranchStatus("*", 0)
    for name in scalar_branches:
        br = chain.GetBranch(name)
        if br is None:
            raise ValueError(
                f"Branch '{name}' not found in tree '{tree_name}'")
        chain.SetBranchStatus(name, 1)
    if array_branch:
        if chain.GetBranch(array_branch) is None:
            raise ValueError(
                f"Branch '{array_branch}' not found in tree '{tree_name}'")
        chain.SetBranchStatus(array_branch, 1)

    # Set up branch addresses for scalars using array.array
    # (numpy buffers don't work with PyROOT for all types)
    buffers = {}
    np_dtypes = {}
    for name in scalar_branches:
        leaf = chain.GetBranch(name).GetLeaf(name)
        type_name = leaf.GetTypeName()

        entry = _TYPE_MAP.get(type_name)
        if entry is None:
            print(f"Warning: skipping branch '{name}' (unsupported type '{type_name}')")
            continue
        typecode, dt = entry
        buf = _array.array(typecode, [0])
        np_dtypes[name] = dt

        chain.SetBranchAddress(name, buf)
        buffers[name] = buf

    # Set up array branch (detect TArrayF vs TArrayS from branch)
    if array_branch:
        br_class = chain.GetBranch(array_branch).GetClassName()
        if "TArrayS" in br_class:
            arr_obj = ROOT.TArrayS()
            arr_dtype = np.int16
        else:
            arr_obj = ROOT.TArrayF()
            arr_dtype = np.float32
        chain.SetBranchAddress(array_branch, arr_obj)
        # Read first entry to determine waveform size
        chain.GetEntry(0)
        wf_size = arr_obj.GetSize()
        waveforms = np.empty((n_to_read, wf_size), dtype=arr_dtype)
    else:
        waveforms = None

    # Pre-allocate scalar output arrays
    scalar_data = {
        name: np.empty(n_to_read, dtype=np_dtypes[name])
        for name in buffers
    }

    read_count = 0
    entry_idx = 0

    while read_count < n_to_read:
        if entry_idx >= n_total:
            break

        nb = chain.GetEntry(entry_idx)
        if nb <= 0:
            entry_idx += 1
            continue

        for name, buf in buffers.items():
            scalar_data[name][read_count] = buf[0]

        if array_branch:
            # Bulk copy via buffer protocol instead of per-element At()
            waveforms[read_count] = np.frombuffer(arr_obj.GetArray(),
                                                  dtype=arr_dtype,
                                                  count=wf_size)

        read_count += 1
        entry_idx += 1

    # Trim arrays if we read fewer than expected
    if read_count < n_to_read:
        for name in scalar_data:
            scalar_data[name] = scalar_data[name][:read_count]
        if waveforms is not None:
            waveforms = waveforms[:read_count]

    features_df = pd.DataFrame(scalar_data)

    if array_branch:
        return features_df, waveforms
    else:
        return features_df
