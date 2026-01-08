{
  description = "ROOT Analysis Development Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    utils.url = "github:ewtodd/Analysis-Utilities";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        analysis-utils = utils.packages.${system}.default;
      in
      {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
            clang-tools
          ];
          buildInputs = with pkgs; [
            analysis-utils
            root
          ];

          shellHook = ''
            echo "ROOT version: $(root-config --version)"
            echo "Analysis-Utilities version: ${analysis-utils.version}"

            STDLIB_PATH="${pkgs.stdenv.cc.cc}/include/c++/${pkgs.stdenv.cc.cc.version}"
            STDLIB_MACHINE_PATH="$STDLIB_PATH/x86_64-unknown-linux-gnu"

            ROOT_INC="$(root-config --incdir)"
            # Local first, then remote, then others
            export CPLUS_INCLUDE_PATH="$PWD/include:$STDLIB_PATH:$STDLIB_MACHINE_PATH:${analysis-utils}/include:$ROOT_INC''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"

            export PKG_CONFIG_PATH="${analysis-utils}/lib/pkgconfig:$PKG_CONFIG_PATH"

            export ROOT_INCLUDE_PATH="$PWD/include:${analysis-utils}/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            # Local lib first means linker will use it preferentially
            export LD_LIBRARY_PATH="$PWD/lib:${analysis-utils}/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };
      }
    );
}
