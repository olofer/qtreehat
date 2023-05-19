# qtreehat
Basic quadtree (2D space) hierarchically arranging quadrupole information for a set of weights and coordinates, with the ability to evaluate all-to-all meshless either $log(r)$ (or $-1/r$) potentials between points (same quadtree, different evaluations). 

## Tests & demos
- `test_qtreehat.m` that verifies the operation of the accuracy parameter
- `viz_qtreehat.m` makes visualizations (impulse and random fields)
- `sim_qtreehat.m` runs a leapfrog integrated $N$-body simulation

## Tasks
- `qtreehat.c` interface to optionally take target points as arguments
- separate tests of the value, gradient pair via random line integration
- split up source code such that the MEX file is a special case application
- use within a $N$-body WASM/C++/JS interactive animation
