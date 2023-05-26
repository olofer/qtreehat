# qtreehat
Basic quadtree (2D space) hierarchically arranging quadrupole information for a set of weights and coordinates. The `MEX` interface has the ability to evaluate either all-to-all or all-to-targets meshless summations, for either $log(r)$ (or $-1/r$) potentials between points (same quadtree, different evaluations). The calculation of the field gradients is always done as a side result.

## Tests & demos
- `test_qtreehat.m` that verifies the operation of the accuracy parameter
- `viz_qtreehat.m` makes visualizations (impulse and random fields)
- `sim_qtreehat.m` runs a leapfrog integrated $N$-body simulation (energy conservation test)
- `line_qtreehat.m` tests/shows potential delta calc. via gradient line integral
- `sample_nnodes.m` counts/shows the number of quadtree nodes required 

## Tasks
- make a non-`MEX` smoke test executable (built via `gcc -Wall`); test nearest-neighbor indexing
- verify separate test executable with `valgrind` 
- use within a $N$-body WASM/C++/JS interactive animation
