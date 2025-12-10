{
  description = "ROOT Analysis Development Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    utils.url = "github:ewtodd/Nuclear-Measurement-Utilities";
    utils-local.url = "path:/path/to/your/local/clone";
  };

  outputs = { self, nixpkgs, flake-utils, utils, utils-local }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        nm-utils-remote = utils.packages.${system}.default;
        nm-utils-local = utils-local.packages.${system}.default;
      in {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [ pkg-config gnumake clang-tools ];
          buildInputs = with pkgs; [ nm-utils root ];

          shellHook = ''
            echo "ROOT version: $(root-config --version)"
            echo "Remote NM Utils: ${nm-utils-remote}"
            echo "Local NM Utils: ${nm-utils-local}"

            STDLIB_PATH="${pkgs.stdenv.cc.cc}/include/c++/${pkgs.stdenv.cc.cc.version}"
            STDLIB_MACHINE_PATH="$STDLIB_PATH/x86_64-unknown-linux-gnu"

            ROOT_INC="$(root-config --incdir)"
            # Local first, then remote, then others
            export CPLUS_INCLUDE_PATH="${nm-utils-local}/include:$PWD/include:$STDLIB_PATH:$STDLIB_MACHINE_PATH:${nm-utils-remote}/include:$ROOT_INC''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"

            export PKG_CONFIG_PATH="${nm-utils-local}/lib/pkgconfig:${nm-utils-remote}/lib/pkgconfig:$PKG_CONFIG_PATH"

            export ROOT_INCLUDE_PATH="${nm-utils-local}/include:$PWD/include:${nm-utils-remote}/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            # Local lib first means linker will use it preferentially
            export LD_LIBRARY_PATH="${nm-utils-local}/lib:$PWD/lib:${nm-utils-remote}/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };
      });
}
