{ pkgs ? import <nixpkgs> { } }:
pkgs.mkShell {
  packages = with pkgs; [ rustup pkg-config fontconfig freetype llvmPackages.llvm python3 pngquant ];
}
