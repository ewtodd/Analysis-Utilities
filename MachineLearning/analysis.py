import pandas as pd
import h5py
import numpy as np
from analyze_psd import (
    regress_waveforms,
    compute_and_analyze_psd,
    process_waveforms,
)
import pickle


def main():

    sources = [0, 2, 3]  # Use integers instead of strings
    merged_file = "all_converted_waveforms.h5"

    print("Loading features...")
    features = pd.read_hdf(merged_file, key="features")

    # Check what source_id values actually exist
    print("Unique source_id values:", features["source_id"].unique())
    print("Source_id data type:", features["source_id"].dtype)

    zero_features = features[features["source_id"] == 0].reset_index(
        drop=True
    )  # Integer 0
    one_features = features[features["source_id"] == 2].reset_index(
        drop=True
    )  # Integer 2

    print(f"Alpha features: {len(zero_features)} events")
    print(f"Gamma/beta features: {len(one_features)} events")

    # Load waveforms from the merged file and filter by indices
    print("Loading waveforms...")
    with h5py.File(merged_file, "r") as f:
        all_waveforms = f["raw_waveform"][:]

    # Get indices for each source
    zero_indices = features[features["source_id"] == 0].index.values
    one_indices = features[features["source_id"] == 2].index.values

    # Extract waveforms for each source
    zero_waveforms_array = all_waveforms[zero_indices]
    one_waveforms_array = all_waveforms[one_indices]

    # Convert to DataFrames
    zero_waveforms = pd.DataFrame(
        zero_waveforms_array,
        columns=[f"sample_{i}" for i in range(zero_waveforms_array.shape[1])],
    )
    one_waveforms = pd.DataFrame(
        one_waveforms_array,
        columns=[f"sample_{i}" for i in range(one_waveforms_array.shape[1])],
    )

    print(f"Alpha waveforms: {zero_waveforms.shape}")
    print(f"Gamma/beta waveforms: {one_waveforms.shape}")
    waveforms = (zero_waveforms, one_waveforms)
    features_tuple = (zero_features, one_features)

    (
        (test_zero_waveforms, test_one_waveforms),
        (test_zero_features, test_one_features),
    ) = regress_waveforms(
        waveforms,
        features_tuple,
        process_func=process_waveforms,
        random_state=42,
        model_file="regressor.pkl",
        output_dir="psd_analysis",
    )

    compute_and_analyze_psd(
        (test_zero_waveforms, test_one_waveforms),
        (test_zero_features, test_one_features),
        output_dir="psd_analysis",
    )


if __name__ == "__main__":
    main()
