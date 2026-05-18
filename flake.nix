{
  description = "Nuclear Measurement Utilities";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
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
        version = "26.5.18";
        toolkit = pkgs.stdenv.mkDerivation {
          pname = "analysis-utilities";
          inherit version;

          src = ./.;

          nativeBuildInputs =
            with pkgs;
            [
              pkg-config
              gnumake
            ]
            ++ pkgs.lib.optionals (!pkgs.stdenv.hostPlatform.isDarwin) [
              autoPatchelfHook
            ];

          buildInputs = with pkgs; [
            root
          ];

          installPhase = ''
            make install PREFIX=$out
          '';

          postFixup = pkgs.lib.optionalString (!pkgs.stdenv.hostPlatform.isDarwin) ''
            for lib in $out/lib/*.so; do
              if [ -f "$lib" ]; then
                patchelf --set-rpath "$out/lib:${pkgs.root}/lib:${pkgs.stdenv.cc.cc.lib}/lib" "$lib" || true
              fi
            done
          '';

          setupHook = ./setup-hook.sh;
        };

        pythonPackage = pkgs.python3Packages.buildPythonPackage {
          pname = "analysis-utilities";
          inherit version;
          src = ./python;
          pyproject = true;

          nativeBuildInputs = [ pkgs.python3Packages.setuptools ];
          postPatch = ''
            substituteInPlace analysis_utilities/__init__.py \
              --replace-fail '@VERSION@' "${version}"
            substituteInPlace pyproject.toml \
              --replace-fail '@VERSION@' "${version}"
          '';
          propagatedBuildInputs = with pkgs.python3Packages; [
            numpy
            pandas
          ];
          doCheck = false;
        };
      in
      {
        packages.default = toolkit;
        packages.pythonPackage = pythonPackage;

        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            root
            gnumake
            pkg-config
            clang-tools
            (python3.withPackages (
              python-pkgs: with python-pkgs; [
                numpy
                pandas
                pytest
                pythonPackage
              ]
            ))
          ];

          shellHook = ''
            export SHELL="${pkgs.bash}/bin/bash"
            echo "Development environment for working on the analysis utilities source"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${pkgs.root}/include"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
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
        python = {
          path = ./templates/python;
          description = "Python analysis development environment with analysis-utils and ML libraries.";
          welcomeText = ''
            Run `nix develop` to enter the development environment.
            The analysis-utils Python package and PlottingUtils bridge are available out of the box.
          '';
        };
      };
    };
}
