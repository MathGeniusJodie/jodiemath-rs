{ pkgs ? import <nixpkgs> { } }:
pkgs.mkShell {
  packages = with pkgs; [ rustup pkg-config fontconfig freetype llvmPackages.llvm python3 pngquant ];
  # .cargo/config.toml builds fontconfig-sys in dlopen mode (so plain `cargo
  # build` needs no pkg-config); the plot-writing lib tests load it from here.
  LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [ pkgs.fontconfig ];
}
