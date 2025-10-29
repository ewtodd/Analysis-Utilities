{
  description = "Digital Signal Processing Framework";

  inputs = { nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable"; };

  outputs = { self, nixpkgs, flake-utils }:
    let
      pkgs = nixpkgs.legacyPackages.${system};

      toolkit = pkgs.stdenv.mkDerivation {
        pname = "dsp-toolkit";
        version = "0.1";

        src = ./.;

        nativeBuildInputs = with pkgs; [ pkg-config ];

        buildInputs = with pkgs; [ root gsl ];

        installPhase = ''
          mkdir -p $out/{lib,include,share/macros}

          cp lib/*.so $out/lib/ 2>/dev/null || true
          cp lib/*.a $out/lib/ 2>/dev/null || true

          if [ -d include ] && [ -n "$(ls -A include/*.h 2>/dev/null)" ]; then
            cp include/*.h $out/include/
          fi

          if [ -d macros ] && [ -n "$(ls -A macros/*.C 2>/dev/null)" ]; then
            cp macros/*.C $out/share/macros/
          fi
        '';
        shellHook = ''
          echo "Digital Signal Processing Development Environment for CAEN digitizers"
          echo "ROOT version: $(root-config --version)"
          echo ""

          # For clang?
          export ROOT_INCLUDE_PATH="$PWD/include:$ROOT_INCLUDE_PATH"
          export LD_LIBRARY_PATH="$PWD/build:$LD_LIBRARY_PATH"
        '';
      };

    in {
      packages.default = toolkit;
      templates = {
        default = {
          path = ./templates/default;
          description = "Standard ROOT waveform analysis pipeline.";
          welcomeText = ''
            Run `nix develop` to enter the development environment.
          '';
        };
      };
    };
}
