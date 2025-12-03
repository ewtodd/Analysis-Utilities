{
  description = "ROOT Analysis Development Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    utils.url = "github:ewtodd/Nuclear-Measurement-Utilities";
  };

  outputs = { self, nixpkgs, flake-utils, utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        nm-utils = utils.packages.${system}.default;
      in {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [ pkg-config gnumake clang-tools ];
          buildInputs = with pkgs; [ nm-utils root ];

          shellHook = ''
            echo "ROOT Waveform Analysis Framework"
            echo "ROOT version: $(root-config --version)"
            echo "Nuclear Measurement Utilities: ${nm-utils}"
            echo ""

            STDLIB_PATH="${pkgs.stdenv.cc.cc}/include/c++/${pkgs.stdenv.cc.cc.version}"
            STDLIB_MACHINE_PATH="$STDLIB_PATH/x86_64-unknown-linux-gnu"

            ROOT_INC="$(root-config --incdir)"
            export CPLUS_INCLUDE_PATH="$STDLIB_PATH:$STDLIB_MACHINE_PATH:${nm-utils}/include:$ROOT_INC''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"

            export PKG_CONFIG_PATH="${nm-utils}/lib/pkgconfig:$PKG_CONFIG_PATH"

            export ROOT_INCLUDE_PATH="${nm-utils}/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            export LD_LIBRARY_PATH="${nm-utils}/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

            if pkg-config --exists nm-utils; then
              echo "Nuclear Measurement Toolkit pkg-config found"
              echo "Includes: $(pkg-config --cflags nm-utils)"
              echo "Libs: $(pkg-config --libs nm-utils)"
            fi
          '';
        };
      });
}

