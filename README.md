# qtreehat
Basic quadtree (2D space) hierarchically arranging quadrupole information for a set of weights and coordinates. The `MEX` interface has the ability to evaluate either all-to-all or all-to-targets meshless summations, for either $log(r)$ (or $-1/r$) potentials between points (same quadtree, different evaluations). The calculation of the field gradients is always done as a side result.

## Tests & demos
- `test_qtreehat.m` that verifies the operation of the accuracy parameter
- `viz_qtreehat.m` makes visualizations (impulse and random fields)
- `sim_qtreehat.m` runs a leapfrog integrated $N$-body simulation
- `sample_nnodes.m` counts/shows the number of quadtree nodes required 

## Tasks
- separate tests of the value, gradient pair via random line integration
- split up source code such that the `MEX` file is a special case application
- use within a $N$-body WASM/C++/JS interactive animation
