from processing import pipeline_get_waveforms, get_file_info, load_specific_features

filepaths = ["../macros/"]
output_file = pipeline_get_waveforms(
    filepaths,
    verbose=True,
    apply_cuts=True,  # Only include events that pass cuts
    max_waveforms=None,  # Limit per file for testing
)

# Load and inspect the converted data
info = get_file_info(output_file)
print("File info:", info)

# Load specific features
features = load_specific_features(
    output_file, columns=["pulse_height", "light_output_keVee", "charge_comparison_psd"]
)
