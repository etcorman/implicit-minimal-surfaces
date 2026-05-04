#  Implicit Minimal Surfaces for Bijective Correspondences

Official MATLAB implementation of the paper **"Implicit Minimal Surfaces for Bijective Correspondences"** *SIGGRAPH 2026 (ACM Transactions on Graphics)*

It computes bijective correspondences between two genus 0 surfaces.

## Usage
The algorithm is launched with the script `run_minimal_surfaces.m`. It will load an `A.obj` and `B.obj` from the folder `Mesh/[name]/` and a set of constraints described in `.pinned` files. The output correspondences are exported in `A_inter.obj` and `B_inter.obj` files in the folder `Results/[name]/` with shared triangulations. `A_tri.obj` and `B_tri.obj` also share triangulation but the UV maps contain the triangle indices of the initial meshes.

## Alternative Implementations
- **C++:** Coming soon
- **Pyton:** (https://github.com/RobinMagnet/implicit-minimal-surfaces)
