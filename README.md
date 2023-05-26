# qtreehat
Basic quadtree (2D space) hierarchically arranging quadrupole information for a set of weights and coordinates. The `MEX` interface has the ability to evaluate either all-to-all or all-to-targets meshless summations, for either $log(r)$ (or $-1/r$) potentials between points (same quadtree, different evaluations). The calculation of the field gradients is always done as a side result. 

## Build & sanity check
Run `./build-run-tests.sh` to compile and run basic smoke tests. Only checked under `WSL2` and requires both `gcc` and `clang` to be installed on the system, as well as `octave`.

## Tests & demos
- `test_qtreehat.m` that verifies the operation of the accuracy parameter
- `viz_qtreehat.m` makes visualizations (impulse and random fields)
- `sim_qtreehat.m` runs a leapfrog integrated $N$-body simulation (energy conservation test)
- `line_qtreehat.m` tests/shows potential delta calc. via gradient line integral
- `sample_nnodes.m` counts/shows the number of quadtree nodes required 

## Tasks
- (kernel-independent) FMM extension?
- C++ interface?
- use within a $N$-body WASM/C++/JS interactive animation (SPH + gravity)
- can it be extended for Helmholtz problems?

## References
This treecode is an extension of the quadtree used here:
- https://github.com/olofer/yasph
