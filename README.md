# qtreehat
Basic quadtree (2D space) hierarchically arranging quadrupole information for a set of weights and coordinates. The `MEX` interface has the ability to evaluate either all-to-all or all-to-targets meshless summations, for either $\sim \log(r)$ (or $\sim -1/r$) potentials between points (same quadtree, different evaluations). The calculation of the field gradients is always done as a side result. 

## Build & sanity check
Run `./build-run-test.sh` to compile and run basic smoke tests. Only checked under `WSL2`. Requires both `gcc` and `clang` to be installed on the system, as well as `valgrind` and `octave`. To build only the `MEX` interface from within `octave` (or `matlab`, but this is not tested) type `build_mex('qtreehat');`. 

Run `./build-run-test++.sh` to build and execute smoke tests of the `C++` interface to the tree code. This requires `g++` and `clang++`. Only checked under `WSL2`.

## Functional tests & demos
- `test_qtreehat.m` that verifies the operation of the accuracy parameter
- `viz_qtreehat.m` makes visualizations (impulse and random fields)
- `sim_qtreehat.m` runs a leapfrog integrated $N$-body simulation (energy conservation test)
- `line_qtreehat.m` tests/shows potential delta calc. via gradient line integral
- `sample_nnodes.m` counts/shows the number of quadtree nodes required 
- `pfgm_qtreehat.m` crude demonstration of Poisson Flow Generative Modeling (in 1D)

## Tasks
- `WIP`: $N$-body WASM/C++/JS interactive animation (SPH + gravity)
- (kernel-independent) FMM extension?
- Octtree variant of the code (at some future time, maybe)

## References
This treecode is an extension of the quadtree used here:
- https://github.com/olofer/yasph
