{
  # TODO: Add a description
  description = "Description of project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
  };
  nixConfig = {
    extra-trusted-public-keys =
      "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, devenv, systems, ... }@inputs:
    let forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in {
      packages = forEachSystem (system: {
        devenv-up = self.devShells.${system}.default.config.procfileScript;
      });

      devShells = forEachSystem (system:
        let pkgs = nixpkgs.legacyPackages.${system};
        in {
          default = devenv.lib.mkShell {
            inherit inputs pkgs;
            modules = [{
              packages = with pkgs; [
                stdenv.cc.cc.lib # required by Jupyter
                zlib
                glibc
                python311Packages.pip
                ruff
                nodePackages.pyright
              ];

              env = { GREET = "󱄅 Nix"; };

              scripts.hello.exec = "echo $GREET";

              enterShell = ''
                hello
              '';

              languages.python = {
                package = pkgs.python311;
                enable = true;
                poetry = {
                  enable = true;
                  activate.enable = true;
                  install.enable = true;
                  install.installRootPackage = true;
                  install.allExtras = true;
                };
              };

              pre-commit.hooks = {
                ruff.enable = true;
                pyright.enable = true;
                shellcheck.enable = true;
                markdownlint.enable = true;
                alejandra.enable = true;
                editorconfig-checker.enable = true;
              };

              dotenv.enable = true;
            }];
          };
        });
    };
}
