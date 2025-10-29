{
  description = "ROOT Waveform Analysis Framework";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    toolkit.url = "github:ewtodd/CAEN-DSP-Toolkit";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system};

      in {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            toolkit.packages.${system}.default
            root
            gsl
            gnumake
            pkg-config
            clang-tools
            (python3.withPackages (python-pkgs:
              with python-pkgs; [
                matplotlib
                numpy
                pandas
                tables
                seaborn
                scikit-learn
                mplhep
                awkward
                pynvim
                ipykernel
                cairosvg
                plotly
                kaleido
                pyarrow
                uproot
                h5py
              ]))
          ];

          shellHook = ''
            echo "ROOT Waveform Analysis Development Environment"
            echo "ROOT version: $(root-config --version)"
            echo ""
            echo "Available commands:"
            echo "  make                 - Build the project"
            echo "  make run-initial     - Run initial waveform processing"
            echo "  make run-calibrate   - Run calibration and plotting"
            echo "  make run-background  - Background subtraction"
            echo "  make run-plots       - Plotting calibrated histograms and average waveforms"
            echo "  make run-PSD         - CC and SI PSD histograms and plotting"
            echo "  make run-optimize    - CC gate optimization"

            # Set up environment
            export ROOT_INCLUDE_PATH="$PWD/include:$ROOT_INCLUDE_PATH"
            export LD_LIBRARY_PATH="$PWD/build:$LD_LIBRARY_PATH"
          '';
        };
      });
}
