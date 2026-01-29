{
  description = "Nuclear Measurement Utilities";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    (flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        toolkit = pkgs.stdenv.mkDerivation {
          pname = "analysis-utilities";
          version = "29.01.2026";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            pkg-config
            autoPatchelfHook
            gnumake
          ];

          buildInputs = with pkgs; [
            root
          ];

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
            cat > $out/lib/pkgconfig/analysis-utilities.pc <<EOF
            prefix=$out
            exec_prefix=\''${prefix}
            libdir=\''${exec_prefix}/lib
            includedir=\''${prefix}/include

            Name: analysis-utilities 
            Description: Analysis Utilities for Nuclear Measurements 
            Libs: -L\''${libdir} -lanalysis-utilities
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

          propagatedBuildInputs = [ pkgs.root ];
        };
      in
      {
        packages.default = toolkit;

        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            root
            gnumake
            pkg-config
            clang-tools
          ];

          shellHook = ''
            export SHELL="/run/current-system/sw/bin/bash"
            echo "Development environment for working on the analysis utilities source"

            STDLIB_PATH="${pkgs.stdenv.cc.cc}/include/c++/${pkgs.stdenv.cc.cc.version}"
            STDLIB_MACHINE_PATH="$STDLIB_PATH/x86_64-unknown-linux-gnu"

            # Build include path in correct order: stdlib -> project -> ROOT
            export CPLUS_INCLUDE_PATH="$STDLIB_PATH:$STDLIB_MACHINE_PATH:$PWD/include:$(root-config --incdir):$CPLUS_INCLUDE_PATH"
            export ROOT_INCLUDE_PATH="$PWD/include:$(root-config --incdir)"
            export LD_LIBRARY_PATH="$PWD/lib:$LD_LIBRARY_PATH"
          '';
        };
      }
    ))
    // {
      templates = {
        default = {
          path = ./templates/standard;
          description = "Standard analysis development environment.";
          welcomeText = ''
            Run `nix develop` to enter the development environment.
            If you have local libraries in include/src, use the included Makefile, and run your macros with root -l macro.cpp+.
          '';
        };
        standard = self.templates.default;
      };
    };
}
