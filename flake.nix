{
  description = "Nuclear Measurement Utilities";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";

    # ROOT master, pinned by commit, for the CUDA build only. It has the CUDA
    # evaluator fixes (root-project/root PRs #23263 and #23753) that no 6.40
    # release carries. Bump it deliberately, by editing the commit here, and run
    # `nix run .#cuda-regression` afterwards.
    root-src = {
      url = "github:root-project/root/e0e5723692467ac53dcc5871e990a68f1e5adf30";
      flake = false;
    };
    # The clad tag that ROOT commit builds against
    # (interpreter/cling/tools/plugins/clad/CMakeLists.txt).
    clad-src = {
      url = "github:vgvassilev/clad/v2.5";
      flake = false;
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      root-src,
      clad-src,
    }:
    (flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaCapabilities = [ "8.9" ];
            cudaForwardCompat = false;
          };
        };

        version = "26.10.22";

        # ROOT for both builds: nixpkgs' 6.40 derivation rebuilt from the
        # root-src master commit above. Released ROOT up to 6.40.04 computes a
        # wrong likelihood for ranged RooSimultaneous fits (root-project/root
        # #23444, fixed in master by 96d8760a77), so no package here uses it.
        rootMaster = pkgs.root.overrideAttrs (old: {
          version = "6.41.01-unstable-2026-10-19";
          src = root-src;

          # Master needs LLVM 22, newer than the patched Clang nixpkgs builds for
          # 6.40, so build ROOT's bundled LLVM and Clang instead and drop
          # nixpkgs' LLVM 20 ones. A derivation-valued attribute is a build input
          # even when unused, hence clang = null.
          clang = null;
          clad_src = clad-src;
          buildInputs = pkgs.lib.filter (
            input: input != null && input != old.clang && input != pkgs.llvm_20
          ) (old.buildInputs or [ ]);

          # nixpkgs deletes every builtins/ directory except two; master keeps
          # mandatory components there, so leave them all. Optional builtins stay
          # off by default, and fail-on-missing still refuses a silent fallback.
          #
          # Master patches clad inside its ExternalProject with `git checkout` and
          # `git apply`, which needs a writable git clone; the clad input is a
          # read-only store path. So apply the same patch to a writable copy here,
          # take it out of ROOT's patch list, and build clad from the copy.
          # --replace-fail stops the build if a ROOT bump changes that list.
          preConfigure = ''
            substituteInPlace cmake/modules/SearchInstalledSoftware.cmake \
              --replace-fail 'set(lcgpackages ' '#set(lcgpackages '
            patchShebangs cmake/unix/

            cladPatch="$PWD/interpreter/cling/tools/plugins/clad/patches/lld_lookup_no_default_path.patch"
            cp -r --no-preserve=mode ${clad-src} "$NIX_BUILD_TOP/clad"
            (cd "$NIX_BUILD_TOP/clad" && git apply --ignore-space-change \
              --ignore-whitespace "$cladPatch")
            substituteInPlace interpreter/cling/tools/plugins/clad/CMakeLists.txt \
              --replace-fail 'list(APPEND _clad_patches_list "lld_lookup_no_default_path.patch")' ""
            appendToVar cmakeFlags "-DCLAD_SOURCE_DIR=$NIX_BUILD_TOP/clad"
          '';

          cmakeFlags = [
            "-Dbuiltin_clang=ON"
            "-Dbuiltin_llvm=ON"
            # Master downloads the HTTP server's civetweb and the web GUI's
            # openui5 at build time, which the sandbox forbids. Nothing here
            # uses either; turning http off also turns webgui off.
            "-Dhttp=OFF"
            "-Dwebgui=OFF"
            # The bundled Clang only searches the FHS include paths, so it cannot
            # find glibc's headers (rootcling fails on <features.h>). nixpkgs'
            # patched Clang for 6.40 is built with the same setting.
            "-DC_INCLUDE_DIRS=${pkgs.lib.getDev pkgs.stdenv.cc.libc}/include"
            "-Dfail-on-missing=ON"
            "-Dfftw3=ON"
            "-Dfitsio=OFF"
            "-Dmathmore=ON"
            "-Dsqlite=OFF"
            "-Dvdt=OFF"
          ];

          # nixpkgs' separate debug output adds -ggdb to every compile, which
          # with the bundled LLVM and Clang means building all of them with
          # debug information.
          separateDebugInfo = false;
        });

        rootWithCuda = rootMaster.overrideAttrs (old: {
          cmakeFlags = old.cmakeFlags ++ [
            "-Dcuda=ON"
            "-DCMAKE_CUDA_ARCHITECTURES=89"
          ];
          buildInputs =
            old.buildInputs
            ++ (with pkgs.cudaPackages; [
              cuda_nvcc
              cuda_cudart
              cccl
            ]);
        });

        # The CUDA build compiles against two of ROOT's private batch-compute
        # headers (RooBatchCompute::Config and CudaStream), which no ROOT
        # install ships. Copying them from the very source rootWithCuda is built
        # from means they can never disagree with its libraries; cp fails loudly
        # if a ROOT bump moves either file.
        rootBatchComputeHeaders =
          pkgs.runCommand "root-batchcompute-headers-${rootWithCuda.version}" { }
            ''
              mkdir -p $out/include/roofit_batchcompute
              cp ${rootWithCuda.src}/roofit/batchcompute/res/RooBatchCompute.h \
                ${rootWithCuda.src}/roofit/batchcompute/src/CudaInterface.h \
                $out/include/roofit_batchcompute/
            '';

        mkUtils =
          { cuda }:
          let
            root = if cuda then rootWithCuda else rootMaster;
            cudaBuildInputs = pkgs.lib.optionals cuda (
              with pkgs.cudaPackages;
              [
                cuda_nvcc
                cuda_cudart
                cccl
              ]
            );
          in
          pkgs.stdenv.mkDerivation {
            pname = "analysis-utilities" + pkgs.lib.optionalString cuda "-cuda";
            inherit version;

            src = ./.;

            nativeBuildInputs =
              with pkgs;
              [
                pkg-config
                cmake
              ]
              ++ pkgs.lib.optionals (!pkgs.stdenv.hostPlatform.isDarwin) [
                autoPatchelfHook
              ];

            buildInputs = [ root ] ++ cudaBuildInputs;

            cmakeFlags = pkgs.lib.optionals cuda [
              "-DAU_USE_CUDA=ON"
              "-DCMAKE_CUDA_ARCHITECTURES=89"
              "-DAU_ROOT_BATCHCOMPUTE_INCLUDE_DIR=${rootBatchComputeHeaders}/include"
            ];

            postFixup = pkgs.lib.optionalString (!pkgs.stdenv.hostPlatform.isDarwin) ''
              for lib in $out/lib/*.so; do
                if [ -f "$lib" ]; then
                  patchelf --set-rpath "$out/lib:${root}/lib:${pkgs.stdenv.cc.cc.lib}/lib" "$lib" || true
                fi
              done
            '';

            setupHook = ./setup-hook.sh;
          };

        utils = mkUtils { cuda = false; };
        utilsCuda = mkUtils { cuda = true; };

        # Regression check (tests/CudaRegression.cpp): CPU/CUDA NLL agreement,
        # RooFit's NLL against a directly computed likelihood, reference NLL
        # values, and the CUDA speed-up over the CPU, on toy photopeak models in
        # the configurations the library fits. It needs a GPU, which the build
        # sandbox does not have, so it is an app rather than a check:
        # `nix run .#cuda-regression`.
        cudaRegression = pkgs.stdenv.mkDerivation {
          pname = "analysis-utilities-cuda-regression";
          inherit version;

          src = ./tests;

          nativeBuildInputs = [ pkgs.makeWrapper ];
          buildInputs = [
            rootWithCuda
            utilsCuda
          ];

          buildPhase = ''
            runHook preBuild
            $CXX -O2 -DAU_ROOFIT_BACKEND_CUDA=1 $(root-config --cflags) \
              -I${utilsCuda}/include CudaRegression.cpp -o cuda-regression \
              -L${utilsCuda}/lib -lanalysis-utils \
              $(root-config --libs) -lRooFit -lRooFitCore
            runHook postBuild
          '';

          # libcuda comes from the host driver, which NixOS exposes under
          # /run/opengl-driver/lib; ROOT_INCLUDE_PATH lets cling find the
          # dictionary header of the custom PDFs.
          installPhase = ''
            runHook preInstall
            install -Dm755 cuda-regression $out/bin/cuda-regression
            wrapProgram $out/bin/cuda-regression \
              --prefix LD_LIBRARY_PATH : /run/opengl-driver/lib \
              --prefix ROOT_INCLUDE_PATH : ${utilsCuda}/include
            runHook postInstall
          '';
        };

        # The published API reference, served straight out of the store by the
        # Caddy vhost on nu. Deliberately delegates to `make docs` so the
        # Doxyfile substitution logic lives in exactly one place.
        docs = pkgs.stdenv.mkDerivation {
          pname = "analysis-utilities-docs";
          inherit version;

          src = ./.;

          nativeBuildInputs = with pkgs; [
            doxygen
            graphviz
          ];

          dontConfigure = true;

          buildPhase = ''
            runHook preBuild
            export DOXYGEN_AWESOME_CSS="${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css"
            make docs

            # A warning here means a broken \ref, a malformed doc block, or a
            # parameter that no longer exists. Failing the build is what keeps
            # the published site honest, and gives CI the check for free.
            if [ -s docs/doxygen-warnings.log ]; then
              echo "doxygen emitted warnings:" >&2
              cat docs/doxygen-warnings.log >&2
              exit 1
            fi
            runHook postBuild
          '';

          installPhase = ''
            runHook preInstall
            mkdir -p $out
            cp -r docs/html/. $out/
            runHook postInstall
          '';
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

        # Drop-in .clangd for the gpu/ subdir so clangd can parse .cu files
        # against the same arch + version policy that nvcc will use.
        clangdConfigFile = (pkgs.formats.yaml { }).generate "dot-clangd" {
          CompileFlags.Add = [
            "--cuda-gpu-arch=sm_89"
            "--no-cuda-version-check"
          ];
          Diagnostics.Suppress = [
            "no_member"
            "nested_name_spec_non_tag"
            "typename_nested_not_found"
            "template_instantiate_undefined"
          ];
        };
      in
      {
        # CPU-batched build (default). Downstream flakes that just want the
        # speedup from doEval pull this.
        packages.default = utils;
        # CUDA build — RunFit() and friends call EvalBackend::Cuda(). Downstream
        # flakes that opt in pull utils.packages.${system}.cuda.
        packages.cuda = utilsCuda;
        packages.pythonPackage = pythonPackage;
        # Static site for au.ethanwtodd.com; CI publishes this to Pages.
        packages.docs = docs;
        # The ROOT each library build links against. Downstream code must build
        # against the matching one, not pkgs.root: a different ROOT is a
        # different ABI, and released ROOT has the ranged-fit bug noted above.
        packages.root = rootMaster;
        packages.rootCuda = rootWithCuda;
        packages.cudaRegression = cudaRegression;
        apps.cuda-regression = {
          type = "app";
          program = "${cudaRegression}/bin/cuda-regression";
          meta.description = "Check the CUDA backend against the CPU backend for agreement and speed";
        };

        devShells.cpu = pkgs.mkShell {
          buildInputs = with pkgs; [
            rootMaster
            gnumake
            pkg-config
            clang-tools
            # API docs: `make docs` / `cmake --build build --target docs`.
            doxygen
            graphviz
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
            export DOXYGEN_AWESOME_CSS="${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css"
            echo "Development environment for working on the analysis utilities source (CPU)"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${rootMaster}/include"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };

        devShells.default = pkgs.mkShell {
          buildInputs =
            (with pkgs; [
              gnumake
              pkg-config
              clang-tools
              # API docs: `make docs` / `cmake --build build --target docs`.
              doxygen
              graphviz
              (python3.withPackages (
                python-pkgs: with python-pkgs; [
                  numpy
                  pandas
                  pytest
                  pythonPackage
                ]
              ))
            ])
            ++ [ rootWithCuda ]
            ++ (with pkgs.cudaPackages; [
              cuda_nvcc
              cuda_cudart
              cccl
            ]);

          shellHook = ''
            export SHELL="${pkgs.bash}/bin/bash"
            export DOXYGEN_AWESOME_CSS="${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css"
            echo "Development environment for working on the analysis utilities source (CUDA, AU_ROOFIT_BACKEND_CUDA=1)"
            export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
            export AU_ROOT_BATCHCOMPUTE_INCLUDE_DIR="${rootBatchComputeHeaders}/include"
            export CPLUS_INCLUDE_PATH="$PWD/include:$AU_ROOT_BATCHCOMPUTE_INCLUDE_DIR''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${rootWithCuda}/include"
            export LD_LIBRARY_PATH="$PWD/lib:/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            mkdir -p "$PWD/gpu"
            install -m 644 ${clangdConfigFile} "$PWD/gpu/.clangd"
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
