{
  description = "Nuclear Measurements Analysis Toolkit";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    (flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        toolkit = pkgs.stdenv.mkDerivation {
          pname = "nm-toolkit";
          version = "0.3";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            pkg-config
            autoPatchelfHook
            gnumake
          ];

          buildInputs = with pkgs; [ root tomlplusplus ];

          buildPhase = ''
            make
          '';

          installPhase = ''
            mkdir -p $out/{lib,include}

            if [ -d lib ] && [ -n "$(ls -A lib/*.so 2>/dev/null)" ]; then
              cp lib/*.so $out/lib/
            else
              echo "ERROR: No shared libraries found in lib/"
              exit 1
            fi

            if [ -d lib ] && [ -n "$(ls -A lib/*.a 2>/dev/null)" ]; then
              cp lib/*.a $out/lib/
            fi

            if [ -d include ] && [ -n "$(ls -A include/*.hpp 2>/dev/null)" ]; then
              cp include/*.hpp $out/include/
            else
              echo "ERROR: No headers found in include/"
              exit 1
            fi

            mkdir -p $out/lib/pkgconfig
            cat > $out/lib/pkgconfig/nm-toolkit.pc <<EOF
            prefix=$out
            exec_prefix=\''${prefix}
            libdir=\''${exec_prefix}/lib
            includedir=\''${prefix}/include

            Name: nm-toolkit 
            Description: Nuclear Measurements Analysis Toolkit 
            Version: 0.3
            Libs: -L\''${libdir} -lnm-toolkit
            Cflags: -I\''${includedir}
            EOF
          '';

          postFixup = ''
            for lib in $out/lib/*.so; do
              if [ -f "$lib" ]; then
                patchelf --set-rpath "$out/lib:${pkgs.root}/lib:${pkgs.stdenv.cc.cc.lib}/lib" "$lib" || true
              fi
            done
          '';

          propagatedBuildInputs = [ pkgs.root pkgs.tomlplusplus ];
        };
      in {
        packages.default = toolkit;
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            root
            tomlplusplus
            gnumake
            pkg-config
            clang-tools
          ];

          shellHook = ''
            echo "Development environment for working on the nuclear measurement toolkit source"

            # Set up environment for local development
            export ROOT_INCLUDE_PATH="$PWD/include:$(root-config --incdir)"
            export CPLUS_INCLUDE_PATH="$PWD/include:$(root-config --incdir):$CPLUS_INCLUDE_PATH"

            # For building locally
            export LD_LIBRARY_PATH="$PWD/lib:$LD_LIBRARY_PATH"
          '';
        };
      })) // {
        templates = {
          default = {
            path = ./templates/standard;
            description = "Standard ROOT waveform analysis pipeline.";
            welcomeText = ''
              Run `nix develop` to enter the development environment.
            '';
          };
          standard = self.templates.default;
        };
      };
}
