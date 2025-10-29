{
  description = "Digital Signal Processing Framework";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    (flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        toolkit = pkgs.stdenv.mkDerivation {
          pname = "dsp-toolkit";
          version = "0.1";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            pkg-config
            autoPatchelfHook
            gnumake
          ];

          buildInputs = with pkgs; [ root ];

          buildPhase = ''
            make
          '';

          installPhase = ''
            mkdir -p $out/{lib,include}

            # Copy built libraries
            if [ -d lib ] && [ -n "$(ls -A lib/*.so 2>/dev/null)" ]; then
              cp lib/*.so $out/lib/
            else
              echo "ERROR: No shared libraries found in lib/"
              exit 1
            fi

            if [ -d lib ] && [ -n "$(ls -A lib/*.a 2>/dev/null)" ]; then
              cp lib/*.a $out/lib/
            fi

            if [ -d include ] && [ -n "$(ls -A include/*.h 2>/dev/null)" ]; then
              cp include/*.h $out/include/
            else
              echo "ERROR: No headers found in include/"
              exit 1
            fi

            mkdir -p $out/lib/pkgconfig
            cat > $out/lib/pkgconfig/dsp-toolkit.pc <<EOF
            prefix=$out
            exec_prefix=\''${prefix}
            libdir=\''${exec_prefix}/lib
            includedir=\''${prefix}/include

            Name: dsp-toolkit
            Description: Digital Signal Processing Framework for CAEN digitizers
            Version: 0.1
            Libs: -L\''${libdir} -lDSP-Toolkit
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
      in { packages.default = toolkit; })) // {
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
