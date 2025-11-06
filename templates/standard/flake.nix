{
  description = "ROOT Waveform Analysis Framework";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    toolkit.url = "github:ewtodd/Nuclear-Measurement-Toolkit";
  };

  outputs = { self, nixpkgs, flake-utils, toolkit }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        nm-toolkit = toolkit.packages.${system}.default;
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            nm-toolkit
            root
            gnumake
            pkg-config
            clang-tools
          ];

          shellHook = ''
            echo "ROOT Waveform Analysis Framework"
            echo "ROOT version: $(root-config --version)"
            echo "Nuclear Measurements Toolkit: ${nm-toolkit}"
            echo ""

            # Make pkg-config aware of the toolkit
            export PKG_CONFIG_PATH="${nm-toolkit}/lib/pkgconfig:$PKG_CONFIG_PATH"

            # Add toolkit to ROOT's search paths
            export ROOT_INCLUDE_PATH="${nm-toolkit}/include:$ROOT_INCLUDE_PATH"
            export LD_LIBRARY_PATH="${nm-toolkit}/lib:$LD_LIBRARY_PATH"

            # Verify toolkit is available
            if pkg-config --exists nm-toolkit; then
              echo "Nuclear Measurement Toolkit pkg-config found"
              echo "Includes: $(pkg-config --cflags nm-toolkit)"
              echo "Libs: $(pkg-config --libs nm-toolkit)"
            fi
          '';
        };
      });
}
