"""Generic TTree to numpy/pandas loader for waveform analysis."""

import numpy as np
import pandas as pd
import ROOT


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

    # Set up branch addresses for scalars
    buffers = {}
    for name in scalar_branches:
        br = chain.GetBranch(name)
        if br is None:
            raise ValueError(
                f"Branch '{name}' not found in tree '{tree_name}'")
        leaf = br.GetLeaf(name)
        type_name = leaf.GetTypeName()

        if type_name in ("Float_t", "float"):
            buf = np.zeros(1, dtype=np.float32)
        elif type_name in ("Double_t", "double"):
            buf = np.zeros(1, dtype=np.float64)
        elif type_name in ("Int_t", "int"):
            buf = np.zeros(1, dtype=np.int32)
        elif type_name in ("UInt_t", "unsigned int"):
            buf = np.zeros(1, dtype=np.uint32)
        elif type_name in ("Short_t", "short"):
            buf = np.zeros(1, dtype=np.int16)
        elif type_name in ("Bool_t", "bool"):
            buf = np.zeros(1, dtype=np.bool_)
        else:
            buf = np.zeros(1, dtype=np.float64)

        chain.SetBranchAddress(name, buf)
        buffers[name] = buf

    # Set up array branch
    if array_branch:
        arr_obj = ROOT.TArrayF()
        chain.SetBranchAddress(array_branch, arr_obj)

    # Read data
    scalar_data = {
        name: np.empty(n_to_read, dtype=buf.dtype)
        for name, buf in buffers.items()
    }
    waveform_list = [] if array_branch else None

    read_count = 0
    entry_idx = 0

    while read_count < n_to_read:
        real_entry = entry_idx

        if real_entry >= n_total:
            break

        nb = chain.GetEntry(real_entry)
        if nb <= 0:
            entry_idx += 1
            continue

        for name, buf in buffers.items():
            scalar_data[name][read_count] = buf[0]

        if array_branch:
            size = arr_obj.GetSize()
            waveform = np.array([arr_obj.At(i) for i in range(size)],
                                dtype=np.float32)
            waveform_list.append(waveform)

        read_count += 1
        entry_idx += 1

    # Trim arrays if we read fewer than expected
    if read_count < n_to_read:
        for name in scalar_data:
            scalar_data[name] = scalar_data[name][:read_count]

    features_df = pd.DataFrame(scalar_data)

    if array_branch:
        waveforms = np.array(waveform_list, dtype=np.float32)
        return features_df, waveforms
    else:
        return features_df
