import os
import ROOT
import h5py
import numpy as np
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from joblib import parallel_backend, Parallel, delayed
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

from PyPlottingUtils import PyPlottingUtils  # Import your class

plot_utils = PyPlottingUtils()
ROOT.gROOT.SetBatch(True)

N_JOBS = 32


def normalize_waveform(waveform):
    """Normalize a waveform by dividing by its maximum value (unless max is zero)."""
    max_val = np.max(waveform)
    if max_val != 0:
        return waveform / max_val
    else:
        return waveform


def process_waveforms(waveform_df, n_jobs=N_JOBS):
    """
    Normalizes each waveform by its maximum value.
    Input:
      waveform_df: DataFrame with shape (n_waveforms, n_samples)
    Output:
      normalized_waveforms: numpy array with shape (n_waveforms, n_samples)
    """
    # Convert DataFrame to numpy array for processing
    waveform_array = waveform_df.values

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        normalized = np.array(list(executor.map(normalize_waveform, waveform_array)))
    return normalized


def regress_waveforms(
    waveforms,
    features,
    process_func,
    random_state=42,
    model_file="regressor.pkl",
    output_dir="psd_analysis",
):
    os.makedirs(output_dir, exist_ok=True)

    energy_lower_dict = {
        "0": 500,
        "1": 0,
    }
    energy_upper_dict = {
        "0": 1750,
        "1": 1750,
    }

    zero_waveforms, one_waveforms = waveforms
    zero_features, one_features = features

    # Filter by energy for training ONLY
    energymask_zero = (
        zero_features["light_output_keVee"] <= energy_upper_dict["0"]
    ) & (zero_features["light_output_keVee"] >= energy_lower_dict["0"])

    energymask_one_beta = (
        one_features["light_output_keVee"] <= energy_upper_dict["1"]
    ) & (one_features["light_output_keVee"] >= energy_lower_dict["1"])

    print(zero_features[energymask_zero])
    # Apply energy masks to DataFrames using .loc
    zero_masked_waveforms = zero_waveforms.loc[energymask_zero].reset_index(drop=True)
    one_masked_waveforms = one_waveforms.loc[energymask_one_beta].reset_index(drop=True)

    frac = 0.30
    n_zero_train = int(len(zero_masked_waveforms) * frac)
    n_one_train = int(len(one_masked_waveforms) * frac)
    min_samples = min(n_zero_train, n_one_train)
    max_samples = 10000
    if min_samples > max_samples:
        print(f"Min samples exceeds maximum desired number. Using {max_samples}...")
        min_samples = max_samples

    print(f"{frac*100}% of zero masked: {n_zero_train}")
    print(f"{frac*100}% of one/beta masked: {n_one_train}")
    print(f"Samples for balanced training: {min_samples}")

    # Sample training data directly from DataFrames
    np.random.seed(random_state)
    train_zero_waveforms = zero_masked_waveforms.sample(
        n=min_samples, random_state=random_state
    )
    train_one_waveforms = one_masked_waveforms.sample(
        n=min_samples, random_state=random_state
    )

    # Get the original indices of training samples for later exclusion
    zero_train_original_indices = (
        zero_waveforms.loc[energymask_zero].iloc[train_zero_waveforms.index].index
    )
    one_train_original_indices = (
        one_waveforms.loc[energymask_one_beta].iloc[train_one_waveforms.index].index
    )

    # Process training waveforms
    with Parallel(n_jobs=2) as parallel:
        train_results = parallel(
            delayed(process_func)(group)
            for group in [train_zero_waveforms, train_one_waveforms]
        )

    x_train = np.vstack(train_results)
    y_train = np.array([0] * len(train_results[0]) + [1] * len(train_results[1]))

    # Train or load model
    if os.path.exists(model_file):
        print("Loading existing model.")
        with open(model_file, "rb") as file:
            regressor = pickle.load(file)
    else:
        print("Training new model...")
        regressor = RandomForestRegressor(
            n_estimators=250,
            max_depth=30,
            random_state=random_state,
            max_samples=0.632,
            max_features="sqrt",
            n_jobs=-1,
            verbose=1,
        )
        regressor.fit(x_train, y_train)
        with open(model_file, "wb") as file:
            pickle.dump(regressor, file)

    # Get training predictions for histogram
    y_train_pred = regressor.predict(x_train)

    # Create test data by dropping training samples
    test_zero_waveforms = zero_waveforms.drop(zero_train_original_indices).reset_index(
        drop=True
    )
    test_one_waveforms = one_waveforms.drop(one_train_original_indices).reset_index(
        drop=True
    )

    test_zero_features = zero_features.drop(zero_train_original_indices).reset_index(
        drop=True
    )
    test_one_features = one_features.drop(one_train_original_indices).reset_index(
        drop=True
    )

    # Process test waveforms for prediction
    with Parallel(n_jobs=2) as parallel:
        test_results = parallel(
            delayed(process_func)(group)
            for group in [test_zero_waveforms, test_one_waveforms]
        )

    X_test = np.vstack(test_results)
    y_test_pred = regressor.predict(X_test)
    # Plot test scores using ROOT - REPLACE MATPLOTLIB SECTION
    zero_test_pred = y_test_pred[: len(test_zero_waveforms)]
    one_test_pred = y_test_pred[len(test_zero_waveforms) :]

    plot_score_histogram(
        zero_test_pred,
        one_test_pred,
        "Test Set Scores",
        os.path.join(output_dir, "test_score_histogram.pdf"),
    )
    # Add regressor output to features
    test_zero_features["Regressor_Output"] = y_test_pred[: len(test_zero_waveforms)]
    test_one_features["Regressor_Output"] = y_test_pred[len(test_zero_waveforms) :]
    # Add regressor output to features
    test_zero_features["Regressor_Output"] = y_test_pred[: len(test_zero_waveforms)]
    test_one_features["Regressor_Output"] = y_test_pred[len(test_zero_waveforms) :]
    plot_feature_importance_waveform_with_average(
        regressor,  # Your trained model
        (test_zero_waveforms, test_one_waveforms),
        (test_zero_features, test_one_features),
        ["0", "2"],
        output_dir=output_dir,
    )
    return (
        (test_zero_waveforms, test_one_waveforms),
        (test_zero_features, test_one_features),
    )


def compute_and_analyze_psd(test_waveforms, test_features, output_dir):
    """Simplified PSD analysis using pre-computed values for ROC comparison only"""

    test_zero_waveforms, test_one_waveforms = test_waveforms
    test_zero_features, test_one_features = test_features

    # Use pre-computed PSD values - just rename for consistency
    test_zero_features = test_zero_features.copy()
    test_one_features = test_one_features.copy()
    test_zero_features["charge_comparison_psd"] = test_zero_features[
        "charge_comparison_psd"
    ]
    test_one_features["charge_comparison_psd"] = test_one_features[
        "charge_comparison_psd"
    ]
    test_zero_features["si_psd"] = test_zero_features["si_psd"]
    test_one_features["si_psd"] = test_one_features["si_psd"]

    # Only do ROC analysis comparing all three methods
    analyze_all_methods(test_zero_features, test_one_features, output_dir=output_dir)


def verify_waveform_feature_alignment(waveforms_df, features_df, n_check=5):
    """Verify that waveforms and features are still aligned"""
    print("Checking waveform-feature alignment...")

    for i in range(min(n_check, len(waveforms_df))):
        # Get waveform from DataFrame
        waveform_row = waveforms_df.iloc[i]
        waveform = waveform_row.values

        # Get corresponding feature row
        feature_row = features_df.iloc[i]

        # Check if they make sense together
        wf_peak = np.max(waveform)
        stored_peak = feature_row.get("pulse_height", "N/A")

        print(f"Row {i}:")
        print(f"  Waveform peak: {wf_peak:.3f}")
        print(f"  Stored peak: {stored_peak}")
        print(f"  Match: {'✓' if abs(wf_peak - stored_peak) < 0.001 else '✗'}")


def analyze_all_methods(test_zero_features, test_one_features, output_dir):
    """ROC analysis comparing ML, Charge Comparison, and si_psd PSD"""

    # Filter zeros (500-1750 keVee)
    zero_mask = (test_zero_features["light_output_keVee"] >= 500) & (
        test_zero_features["light_output_keVee"] <= 1750
    )
    test_zero_features_filtered = test_zero_features[zero_mask].reset_index(drop=True)

    # Create combined dataset
    test_features = pd.concat(
        [test_zero_features_filtered, test_one_features]
    ).reset_index(drop=True)

    y_true = np.array(
        [0] * len(test_zero_features_filtered) + [1] * len(test_one_features)
    )

    # All three methods for comparison
    all_methods = ["Regressor_Output", "charge_comparison_psd", "si_psd"]
    all_method_names = ["Random Forest", "Charge Comparison", "Shape Indicator"]

    # Create ROC curve plot comparing all three
    plot_unified_roc_curves(
        test_features,
        y_true,
        all_methods,
        all_method_names,
        output_dir=output_dir,
    )


def plot_feature_importance_waveform_with_average(
    regressor,
    test_waveforms,
    test_features,
    sources,
    output_dir="psd_analysis",
):
    """Plot the feature importance and average waveform (both normalized) using ROOT."""
    os.makedirs(output_dir, exist_ok=True)

    # Unpack test data
    test_zero_waveforms, test_one_waveforms = test_waveforms
    test_zero_features, test_one_features = test_features

    # Filter zero events for training energy range (500-1750 keVee) to match training
    zero_mask = (test_zero_features["light_output_keVee"] >= 500) & (
        test_zero_features["light_output_keVee"] <= 1750
    )
    zero_waveforms_filtered = test_zero_waveforms[zero_mask]

    # Process the filtered waveforms to match ML input format
    zero_waveforms_processed = process_waveforms(zero_waveforms_filtered)

    # Calculate average normalized waveform
    avg_waveform = np.mean(zero_waveforms_processed, axis=0)

    # Extract feature importances
    importances = regressor.feature_importances_

    # Ensure importances match waveform length
    if len(importances) != len(avg_waveform):
        print(
            f"Warning: Feature importance length ({len(importances)}) != waveform length ({len(avg_waveform)})"
        )
        # Pad or truncate as needed
        if len(importances) < len(avg_waveform):
            zero_padded_importances = np.zeros_like(avg_waveform)
            zero_padded_importances[: len(importances)] = importances
        else:
            zero_padded_importances = importances[: len(avg_waveform)]
    else:
        zero_padded_importances = importances

    # Normalize both to [0, 1] for comparison
    avg_waveform_norm = (
        avg_waveform / np.max(avg_waveform)
        if np.max(avg_waveform) > 0
        else avg_waveform
    )
    importances_norm = (
        zero_padded_importances / np.max(zero_padded_importances)
        if np.max(zero_padded_importances) > 0
        else zero_padded_importances
    )

    # Time axis (2 ns sampling)
    x_values = np.arange(len(avg_waveform)) * 2

    # Create ROOT canvas
    canvas = ROOT.TCanvas(
        "c_waveform_importance", "Waveform and Feature Importance", 1600, 1000
    )
    plot_utils.ConfigureCanvas(canvas)

    # Create TGraphs
    graph_waveform = ROOT.TGraph(
        len(x_values), x_values.astype(np.float64), avg_waveform_norm.astype(np.float64)
    )
    graph_importance = ROOT.TGraph(
        len(x_values), x_values.astype(np.float64), importances_norm.astype(np.float64)
    )

    # Configure waveform graph
    graph_waveform.SetLineColor(ROOT.kBlue + 1)
    graph_waveform.SetLineWidth(3)
    graph_waveform.SetTitle("")
    graph_waveform.GetXaxis().SetTitle("Time [ns]")
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(0, 1.1)

    # Configure importance graph
    graph_importance.SetLineColor(ROOT.kRed + 1)
    graph_importance.SetLineWidth(3)
    graph_importance.SetLineStyle(2)  # Dashed line

    # Draw graphs
    graph_waveform.Draw("AL")
    graph_importance.Draw("L SAME")

    # Add vertical line at trigger (17 * 2 = 34 ns)
    trigger_line = ROOT.TLine(34, 0, 34, 1.1)
    trigger_line.SetLineColor(ROOT.kGreen + 2)
    trigger_line.SetLineStyle(9)
    trigger_line.SetLineWidth(2)
    trigger_line.Draw()

    # Create legend
    leg = plot_utils.CreateLegend(0.55, 0.6, 0.92, 0.85)
    leg.AddEntry(graph_waveform, "Average #zero Waveform", "l")
    leg.AddEntry(graph_importance, "Feature Importance", "l")
    leg.AddEntry(trigger_line, "Trigger (34 ns)", "l")
    leg.Draw()

    # Save the plot
    output_path = os.path.join(
        output_dir, "feature_importance_waveform_with_average.pdf"
    )
    canvas.SaveAs(output_path)
    canvas.Close()

    print(f"Feature importance plot saved to {output_path}")


def analyze_auc_and_feature_importance_by_sample_size(
    zero_waveforms,
    one_waveforms,
    zero_features,
    one_features,
    sample_sizes,
    process_func,
    output_dir="psd_analysis",
    energy_cut=(500, 1750),
):
    """
    Trains models with different training sample sizes, calculates AUC, and plots feature importance.
    Applies energy filtering to zero events like the original code.

    Args:
        zero_waveforms: DataFrame of zero waveforms
        one_waveforms: DataFrame of one/beta waveforms
        zero_features: DataFrame with zero features
        one_features: DataFrame with one/beta features
        sample_sizes: list of integers with training sample sizes
        process_func: preprocessing function for waveforms
        output_dir: directory to save plots
        energy_cut: tuple (min_keV, max_keV) for zero energy filtering

    Returns:
        auc_list: list of AUC values for each sample size
        feature_importances: list of arrays of feature importances
    """
    importances_all = []
    auc_list = []

    os.makedirs(output_dir, exist_ok=True)

    # Apply energy cut for zero events (like original code does for training)
    zero_mask = (zero_features["light_output_keVee"] >= energy_cut[0]) & (
        zero_features["light_output_keVee"] <= energy_cut[1]
    )
    zero_waveforms_filtered = zero_waveforms.loc[zero_mask].reset_index(drop=True)
    zero_features_filtered = zero_features.loc[zero_mask].reset_index(drop=True)

    print(
        f"After energy filter ({energy_cut[0]}-{energy_cut[1]} keVee): {len(zero_waveforms_filtered)} zero events"
    )
    print(f"Total one/beta events: {len(one_waveforms)}")

    for n_samples in sample_sizes:
        print(f"Training with {n_samples} samples per class (after energy cut)...")

        # Sample the waveforms and features from filtered zero data
        zero_sample = zero_waveforms_filtered.sample(
            n=min(n_samples, len(zero_waveforms_filtered)), random_state=42
        ).reset_index(drop=True)
        one_sample = one_waveforms.sample(
            n=min(n_samples, len(one_waveforms)), random_state=42
        ).reset_index(drop=True)

        # Prepare training data
        x_train_zero = process_func(zero_sample)
        x_train_one = process_func(one_sample)
        X_train = np.vstack((x_train_zero, x_train_one))
        y_train = np.array([0] * len(x_train_zero) + [1] * len(x_train_one))

        # Train model
        from sklearn.ensemble import RandomForestRegressor

        model = RandomForestRegressor(
            n_estimators=250,
            max_depth=30,
            max_features="sqrt",
            random_state=42,
            n_jobs=-1,
        )
        model.fit(X_train, y_train)

        # Save feature importances
        importances_all.append(model.feature_importances_)

        # Prepare test data: Use remaining filtered data excluding train samples
        zero_test = zero_waveforms_filtered.drop(zero_sample.index).reset_index(
            drop=True
        )
        one_test = one_waveforms.drop(gamma_sample.index).reset_index(drop=True)

        x_test_zero = process_func(zero_test)
        x_test_one = process_func(one_test)
        X_test = np.vstack((x_test_zero, x_test_one))
        y_test_true = np.array([0] * len(x_test_zero) + [1] * len(x_test_one))

        y_test_pred = model.predict(X_test)

        from sklearn.metrics import roc_curve, auc

        fpr, tpr, _ = roc_curve(y_test_true, y_test_pred)
        auc_score = auc(fpr, tpr)
        auc_list.append(auc_score)
        print(f"Sample size: {n_samples}, AUC: {auc_score:.4f}")

    # Plot AUC vs Sample size
    plt.figure(figsize=(30, 20))
    plt.plot(sample_sizes, auc_list, marker="o", linestyle="-", color="blue")
    plt.xlabel("Number of Training Samples per Class")
    plt.ylabel("ROC AUC")
    plt.title(
        f"ROC AUC vs Number of Training Samples\n(Alpha energy: {energy_cut[0]}-{energy_cut[1]} keVee)"
    )
    plt.grid(True)
    plt.xscale("log")
    plt.tight_layout()
    auc_plot_path = os.path.join(output_dir, "auc_vs_samples_filtered.pdf")
    plt.savefig(auc_plot_path, dpi=300)
    plt.close()

    # Plot feature importance comparison
    plt.figure(figsize=(30, 20))
    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_sizes)))
    time_axis = np.arange(len(importances_all[0])) * 2  # Assuming 2 ns sampling

    for i, (imp, n_samples) in enumerate(zip(importances_all, sample_sizes)):
        plt.plot(
            time_axis, imp, label=f"{n_samples} samples", color=colors[i], linewidth=3
        )

    plt.xlabel("Time [ns]")
    plt.ylabel("Feature Importance")
    plt.title(
        f"Feature Importance vs Time for Different Training Sample Sizes\n(Alpha energy: {energy_cut[0]}-{energy_cut[1]} keVee)"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fi_plot_path = os.path.join(
        output_dir, "feature_importance_vs_samples_filtered.pdf"
    )
    plt.savefig(fi_plot_path, dpi=300)
    plt.close()

    print(f"Saved AUC plot to {auc_plot_path}")
    plt.close()

    return auc_list, importances_all


def plot_score_histogram(zero_scores, one_scores, title, output_path):
    """Plot score histogram using ROOT"""
    canvas = ROOT.TCanvas("c_scores", title, 1200, 900)
    plot_utils.ConfigureCanvas(canvas, logy=True)

    # Determine range
    all_scores = np.concatenate([zero_scores, one_scores])
    score_min = np.min(all_scores)
    score_max = np.max(all_scores)

    # Create histograms
    h_zero = ROOT.TH1F("h_zero", "", 75, score_min, score_max)
    h_one = ROOT.TH1F("h_one", "", 75, score_min, score_max)

    # Fill histograms
    for val in zero_scores:
        h_zero.Fill(val)
    for val in one_scores:
        h_one.Fill(val)

    # Configure
    plot_utils.ConfigureHistogram(h_zero, ROOT.kRed + 1)
    plot_utils.ConfigureHistogram(h_one, ROOT.kGreen + 2)

    h_zero.GetXaxis().SetTitle("Regressor Output")
    h_zero.GetYaxis().SetTitle("Counts")
    h_zero.SetTitle("")

    # Draw
    max_val = max(h_zero.GetMaximum(), h_one.GetMaximum())
    h_zero.SetMaximum(max_val * 1.2)
    h_zero.Draw("HIST")
    h_one.Draw("HIST SAME")

    # Legend
    leg = plot_utils.CreateLegend(0.65, 0.7, 0.92, 0.85)
    leg.AddEntry(h_zero, f"Am-241 (#zero)", "f")
    leg.AddEntry(h_one, f"Na-22 (#one)", "f")
    leg.Draw()

    canvas.SaveAs(output_path)
    canvas.Close()


def plot_unified_roc_curves(
    test_features, y_true, methods, method_names, output_dir="psd_analysis"
):
    """Plot ROC curves for all methods using ROOT"""
    os.makedirs(output_dir, exist_ok=True)
    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2]
    target_fpr = 0.05

    canvas = ROOT.TCanvas("c_unified_roc", "Unified ROC", 1200, 900)
    plot_utils.ConfigureCanvas(canvas)

    roc_graphs = []
    # Center the legend at the top - adjust these coordinates
    leg = plot_utils.CreateLegend(0.37, 0.2, 0.92, 0.4)
    leg.SetMargin(0.1)

    for i, (method, name) in enumerate(zip(methods, method_names)):
        scores = test_features[method].values

        # Fix for si_psd: invert scores since higher values correspond to negative class
        if method == "si_psd":
            scores_to_use = -scores  # Invert the scores
        else:
            scores_to_use = scores

        fpr, tpr, thresholds = roc_curve(y_true, scores_to_use)

        # Find threshold and TPR at 5% FPR
        index = np.argmin(np.abs(fpr - target_fpr))
        threshold_at_5pct_fpr = thresholds[index]
        tpr_at_5pct_fpr = tpr[index]
        actual_fpr = fpr[index]

        print(f"{name}:")
        if method == "si_psd":
            # For si_psd, show both the inverted threshold and original threshold
            original_threshold = (
                -threshold_at_5pct_fpr
            )  # Convert back to original scale
            print(
                f"  Inverted threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}"
            )
            print(f"  Original threshold (lower is better): {original_threshold:.6f}")
        else:
            print(f"  Threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}")
        print(f"  TPR at {actual_fpr:.3f} FPR: {tpr_at_5pct_fpr:.3f}")

        auc_score = auc(fpr, tpr)

        # Create TGraph
        roc_graph = ROOT.TGraph(len(fpr), fpr, tpr)
        roc_graph.SetLineColor(colors[i])
        roc_graph.SetLineWidth(3)
        roc_graphs.append(roc_graph)

        # Draw
        if i == 0:
            roc_graph.SetTitle("")
            roc_graph.GetXaxis().SetTitle("False Positive Rate (1 - Specificity)")
            roc_graph.GetYaxis().SetTitle("True Positive Rate (Sensitivity)")
            roc_graph.GetXaxis().SetRangeUser(0, 1)
            roc_graph.GetYaxis().SetRangeUser(0, 1)
            roc_graph.Draw("AL")
        else:
            roc_graph.Draw("L SAME")

        leg.AddEntry(roc_graph, f"{name} (AUC = {auc_score:.2f})", "l")

    # Add diagonal line
    diagonal = ROOT.TLine(0, 0, 1, 1)
    diagonal.SetLineColor(ROOT.kBlack)
    diagonal.SetLineStyle(2)
    diagonal.Draw()

    leg.Draw()
    canvas.SaveAs(f"{output_dir}/unified_roc_curves.pdf")
    canvas.SaveAs(f"{output_dir}/unified_roc_curves.png")
