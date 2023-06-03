/*
 * WebAssembly simulation of self-gravitating gas in 2D.
 * Uses quadtree for both fixed-radius neighbor search (pressure) and long-range interaction (gravity). 
 * The (ideal) gas is represented with smoothed particle hydrodynamics.
 * The simulation is lossless (storage: internal energy of gas, gravitational configuration, kinetic energy).
 *
 */

#include <emscripten.h>
#include <cmath>
#include <cstring>
#include "qtreehat.hpp"

const double sph_kernel_h = 3.0;
const int NMAX = 10000;
const int LEAF = 8;

struct tParticle {
  double m;
  double x;
  double y;
  double u;
  double rho;
  double vx;
  double vy;
};

static tParticle particle[NMAX];
static tQuadTreeHandle<NMAX, LEAF> qtree;

extern "C" {

EMSCRIPTEN_KEEPALIVE
double getRandomJS() {
  return emscripten_random();
}

EMSCRIPTEN_KEEPALIVE
void initializeUniformly(int n,
                         double m, 
                         double xmin, 
                         double xmax, 
                         double ymin, 
                         double ymax)
{
  const double xrange = xmax - xmin;
  const double yrange = ymax - ymin;
  for (int i = 0; i < n; i++) {
    particle[i].x = xmin + getRandomJS() * xrange;
    particle[i].y = ymin + getRandomJS() * yrange;
    particle[i].m = m;
    particle[i].vx = 0.0;
    particle[i].vy = 0.0;
  }
}

EMSCRIPTEN_KEEPALIVE
double getParticleX(int i) {
  return particle[i].x;
}

EMSCRIPTEN_KEEPALIVE
double getParticleY(int i) {
  return particle[i].y;
}

EMSCRIPTEN_KEEPALIVE
void leapfrogTimestep(int n) {
  return;
}

}
