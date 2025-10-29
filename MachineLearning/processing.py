import numpy as np
import pandas as pd
import glob
import ROOT
import os
import h5py
from multiprocessing import Pool

##############################################
# Utils
##############################################


def load_data(filepaths, filename_pattern="processed_waveforms.root", verbose=False):
    """Load ROOT files with the specified pattern"""
    root_files = []
    for filepath in filepaths:
        matching_files = glob.glob(os.path.join(filepath, filename_pattern))
        if matching_files:
            root_files.append(matching_files[0])
            if verbose:
                print("Found file:", matching_files[0])
        else:
            # Try looking in subdirectories
            matching_files = glob.glob(
                os.path.join(filepath, "**", filename_pattern), recursive=True
            )
            if matching_files:
                root_files.append(matching_files[0])
                if verbose:
                    print("Found file:", matching_files[0])
    return root_files


def _out_h5_path(base_name, output_dir):
    base, _ = os.path.splitext(os.path.basename(base_name))
    out_h5 = os.path.join(output_dir, f"{base}_converted.h5")
    return out_h5


##############################################
# Extraction and H5 Writing for One File
##############################################
def process_single_file_to_h5(args):
    (root_file, filepath, max_waveforms, verbose, output_dir, apply_cuts) = args

    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        if verbose:
            print(f"Error opening file: {root_file}")
        return None

    tree = f.Get("features")
    if not tree:
        if verbose:
            print(f"Error: TTree 'features' not found in file: {root_file}")
        f.Close()
        return None

    num_entries = tree.GetEntries()
    if num_entries == 0:
        if verbose:
            print(f"Skipping empty file: {root_file}")
        f.Close()
        return None

    # Set up branch variables
    pulse_height = np.array([0.0], dtype=np.float32)
    peak_position = np.array([0], dtype=np.int32)
    trigger_position = np.array([0], dtype=np.int32)
    long_integral = np.array([0.0], dtype=np.float32)
    source_id = np.array([0], dtype=np.int32)
    passes_cuts = np.array([False], dtype=bool)
    negative_fraction = np.array([0.0], dtype=np.float32)
    light_output_keVee = np.array([0.0], dtype=np.float32)
    charge_comparison_psd = np.array([0.0], dtype=np.float32)
    si_psd = np.array([0.0], dtype=np.float32)
    samples = ROOT.TArrayS()

    # Set branch addresses
    tree.SetBranchAddress("pulse_height", pulse_height)
    tree.SetBranchAddress("peak_position", peak_position)
    tree.SetBranchAddress("trigger_position", trigger_position)
    tree.SetBranchAddress("long_integral", long_integral)
    tree.SetBranchAddress("source_id", source_id)
    tree.SetBranchAddress("passes_cuts", passes_cuts)
    tree.SetBranchAddress("negative_fraction", negative_fraction)
    tree.SetBranchAddress("light_output_keVee", light_output_keVee)
    tree.SetBranchAddress("charge_comparison_psd", charge_comparison_psd)
    tree.SetBranchAddress("si_psd", si_psd)
    tree.SetBranchAddress("Samples", samples)

    rows = []
    waveforms = []
    entries_processed = 0
    entries_accepted = 0
    entries_cut_rejected = 0

    for entry_index in range(num_entries):
        if max_waveforms is not None and entries_accepted >= max_waveforms:
            break

        if tree.GetEntry(entry_index) <= 0:
            continue

        entries_processed += 1

        # Apply cuts filter if requested
        if apply_cuts and not passes_cuts[0]:
            entries_cut_rejected += 1
            continue

        # Extract waveform
        waveform_length = samples.GetSize()
        if waveform_length == 0:
            continue

        waveform = np.array(
            [samples.At(i) for i in range(waveform_length)], dtype=np.float32
        )

        # Create feature dictionary
        feature_dict = {
            "pulse_height": pulse_height[0],
            "peak_position": peak_position[0],
            "trigger_position": trigger_position[0],
            "long_integral": long_integral[0],
            "source_id": source_id[0],
            "passes_cuts": passes_cuts[0],
            "negative_fraction": negative_fraction[0],
            "light_output_keVee": light_output_keVee[0],
            "charge_comparison_psd": charge_comparison_psd[0],
            "si_psd": si_psd[0],
            "source_file": os.path.basename(filepath),
        }

        rows.append(feature_dict)
        waveforms.append(waveform)
        entries_accepted += 1

    if verbose:
        print(f"=== Processing Summary for {filepath} ===")
        print(f"Total entries processed: {entries_processed}")
        print(f"Valid waveforms saved: {entries_accepted}")
        if apply_cuts:
            print(f"Entries rejected by cuts: {entries_cut_rejected}")
        if entries_processed > 0:
            acceptance_rate = (entries_accepted / entries_processed) * 100
            print(f"Acceptance rate: {acceptance_rate:.1f}%")
        print(f"=" * 50)

    f.Close()

    if len(rows) == 0:
        return None

    # Convert to arrays
    arr_waveforms = np.array(waveforms, dtype=np.float32)
    features = pd.DataFrame(rows)

    # Save to HDF5
    out_h5 = _out_h5_path(filepath, output_dir)
    with h5py.File(out_h5, "w") as hf:
        hf.create_dataset("raw_waveform", data=arr_waveforms, compression="gzip")
    features.to_hdf(out_h5, key="features", mode="a", format="table")

    if verbose:
        print(f"Saved {arr_waveforms.shape[0]} waveforms/features: {out_h5}")

    return (out_h5, arr_waveforms.shape[0])


##############################################
# Parallel Extraction: Each File to Its Own HDF5
##############################################
def process_and_extract_waveforms_parallel_to_h5(
    root_files,
    filepaths,
    output_dir,
    verbose=True,
    max_waveforms=None,
    apply_cuts=False,
):
    """Process ROOT files in parallel and convert to HDF5 format"""
    os.makedirs(output_dir, exist_ok=True)

    file_args = [
        (rf, fp, max_waveforms, verbose, output_dir, apply_cuts)
        for rf, fp in zip(root_files, filepaths)
    ]

    with Pool(processes=os.cpu_count()) as pool:
        outputs = pool.map(process_single_file_to_h5, file_args)

    valid_outputs = [item for item in outputs if item is not None]
    return valid_outputs


##############################################
# Merge all per-file HDF5 into one
##############################################
def merge_all_files_to_one(output_files, out_h5):
    """Merge multiple HDF5 files into a single file"""
    arr_waveforms = []
    all_features = []

    for h5file, n_events in output_files:
        with h5py.File(h5file, "r") as f:
            waveforms = f["raw_waveform"][:]
        arr_waveforms.append(waveforms)

        features = pd.read_hdf(h5file, key="features", mode="r")
        all_features.append(features)

    # Concatenate all data
    arr_waveforms = np.concatenate(arr_waveforms, axis=0)
    features = pd.concat(all_features, ignore_index=True)

    # Save merged data
    with h5py.File(out_h5, "w") as f:
        f.create_dataset("raw_waveform", data=arr_waveforms, compression="gzip")
    features.to_hdf(out_h5, key="features", mode="a", format="table")


##############################################
# Chunked Loader
##############################################
def load_in_chunks_from_h5(h5_file, chunk_size=5000, features_key="features"):
    """Load data from HDF5 file in chunks"""
    features_store = pd.HDFStore(h5_file, "r")

    with h5py.File(h5_file, "r") as f:
        n_rows = f["raw_waveform"].shape[0]

    for start in range(0, n_rows, chunk_size):
        stop = min(start + chunk_size, n_rows)

        with h5py.File(h5_file, "r") as f:
            waveforms = f["raw_waveform"][start:stop]

        meta = features_store.select(features_key, start=start, stop=stop)
        meta = meta.reset_index(drop=True)

        yield waveforms, meta

    features_store.close()


##############################################
# Top-Level Pipeline
##############################################
def pipeline_get_waveforms(
    filepaths,
    verbose=False,
    max_waveforms=None,
    output_dir="processed_files",
    output_file="all_converted_waveforms.h5",
    apply_cuts=False,
    filename_pattern="processed_waveforms.root",
):
    """
    Main pipeline to convert ROOT files to HDF5 format

    Parameters:
    -----------
    filepaths : list
        List of directories containing ROOT files
    verbose : bool
        Print detailed processing information
    max_waveforms : int or None
        Maximum number of waveforms to process per file
    output_dir : str
        Directory to store individual processed files
    output_file : str
        Final merged output file name
    apply_cuts : bool
        Whether to apply the passes_cuts filter
    filename_pattern : str
        Pattern to match ROOT files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Find ROOT files
    root_files = load_data(filepaths, filename_pattern, verbose)
    if not root_files:
        print("No ROOT files found!")
        return None

    # Check for existing processed files
    all_expected_processed = [
        _out_h5_path(fp, output_dir)
        for fp in filepaths
        if any(rf.startswith(fp) for rf in root_files)
    ]
    exist_mask = [os.path.exists(h5file) for h5file in all_expected_processed]

    if all(exist_mask) and len(exist_mask) > 0:
        if verbose:
            print("All processed files found; only merging.")
        file_outputs = []
        for h5file in all_expected_processed:
            try:
                with h5py.File(h5file, "r") as f:
                    n_events = f["raw_waveform"].shape[0]
                file_outputs.append((h5file, n_events))
            except Exception as e:
                if verbose:
                    print(f"Error reading {h5file}: {e}")

        if file_outputs:
            merge_all_files_to_one(file_outputs, output_file)
            if verbose:
                total_events = sum(n for _, n in file_outputs)
                print(
                    f"All processed data merged to {output_file} (events: {total_events})"
                )
            return output_file

    # Process files
    if verbose:
        print(f"Processing {len(root_files)} ROOT files...")

    file_outputs = process_and_extract_waveforms_parallel_to_h5(
        root_files, filepaths, output_dir, verbose, max_waveforms, apply_cuts
    )

    if not file_outputs:
        print("No files were successfully processed!")
        return None

    # Merge all processed files
    merge_all_files_to_one(file_outputs, output_file)

    if verbose:
        total_events = sum(n for _, n in file_outputs)
        print(f"All processed data merged to {output_file} (events: {total_events})")

    return output_file


##############################################
# Convenience function for loading specific columns
##############################################
def load_specific_features(h5_file, columns=None):
    """
    Load specific feature columns from HDF5 file

    Parameters:
    -----------
    h5_file : str
        Path to HDF5 file
    columns : list or None
        List of column names to load. If None, loads all columns.

    Returns:
    --------
    pandas.DataFrame
        DataFrame with requested columns
    """
    if columns is None:
        return pd.read_hdf(h5_file, key="features")
    else:
        return pd.read_hdf(h5_file, key="features", columns=columns)


def get_file_info(h5_file):
    """
    Get information about the HDF5 file contents

    Parameters:
    -----------
    h5_file : str
        Path to HDF5 file

    Returns:
    --------
    dict
        Dictionary with file information
    """
    info = {}

    with h5py.File(h5_file, "r") as f:
        if "raw_waveform" in f:
            waveform_shape = f["raw_waveform"].shape
            info["n_waveforms"] = waveform_shape[0]
            info["waveform_length"] = (
                waveform_shape[1] if len(waveform_shape) > 1 else None
            )

    try:
        features = pd.read_hdf(h5_file, key="features")
        info["n_features"] = len(features)
        info["feature_columns"] = list(features.columns)
        info["feature_dtypes"] = dict(features.dtypes)
    except Exception as e:
        info["feature_error"] = str(e)

    return info
