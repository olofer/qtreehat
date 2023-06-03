/*
 * WebAssembly simulation of self-gravitating gas in 2D.
 * Uses quadtree for both fixed-radius neighbor search (pressure) and long-range interaction (gravity). 
 * The (ideal) gas is represented with smoothed particle hydrodynamics.
 * The simulation is lossless (storage: internal energy of gas, gravitational configuration, kinetic energy).
 *
 */

// TODO: total system energy (document clearly)
// TODO: variate generation in circular region
// TODO: viscosity loss

#include <emscripten.h>
#include <cmath>
#include <cstring>
#include "qtreehat.hpp"

const double gas_gamma = 5.0 / 3.0;

const double nearfield_epksq = 1.0e-6;

const double sph_kernel_h = 5.0;
const double sph_kernel_radius = sph_kernel_h * 2.0;

// Wendland C2 kernel; support radius = 2 * h
void evaluate_kernel_wc2(double x, 
                         double y, 
                         double h, 
                         double* w,
                         double* wx,
                         double* wy)
{
  const double _one_pi = 3.14159265358979323846;
  const double sigma = 7.0 / (64.0 * _one_pi);
  const double h2 = h * h;
  const double N = sigma / h2;
  const double q2 = (x * x + y * y) / h2;
  if (q2 > 4.0) {
    *w = 0.0;
    *wx = 0.0;
    *wy = 0.0;
    return;
  }
  const double q = sqrt(q2);
  const double tmp = 4.0 - 4.0 * q + q2;
  const double W = N * (1.0 + 2.0 * q) * tmp * tmp;
  *w = W;
  const double Wprime = -10.0 * (N / h2) * (2.0 - q) * tmp;
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

const int NMAX = 10000;
const int LEAF = 8;

struct tParticle {
  double m;
  double x;
  double y;
  double u;
  double rho;
  double p;
  double vx;
  double vy;
  double phi;
  double udot;
  double vxdot;
  double vydot;
};

static double time;
static tParticle particle[NMAX];
static tQuadTreeHandle<NMAX, LEAF> qtree;

void density_summation_callback(int i, int j, void* aux) {
  tParticle* p = reinterpret_cast<tParticle*>(aux);
  const double xi = p[i].x;
  const double yi = p[i].y;
  const double xj = p[j].x;
  const double yj = p[j].y;
  double w, wx, wy;
  evaluate_kernel_wc2(xi - xj, yi - yj, sph_kernel_h, &w, &wx, &wy);
  p[i].rho += p[j].m * w;
}

void dot_summation_callback(int i, int j, void* aux) {
  if (i == j) return;
  tParticle* p = reinterpret_cast<tParticle*>(aux);
  const double pi = p[i].p;
  const double rhoi = p[i].rho;
  const double xi = p[i].x;
  const double yi = p[i].y;
  const double pj = p[j].p;
  const double rhoj = p[j].rho;
  const double xj = p[j].x;
  const double yj = p[j].y;
  const double dxij = xi - xj;
  const double dyij = yi - yj;
  double w, wx, wy;
  evaluate_kernel_wc2(dxij, dyij, sph_kernel_h, &w, &wx, &wy);
  const double mj = p[j].m;
  const double Ci = pi / (rhoi * rhoi);
  const double Aij = mj * (Ci + pj / (rhoj * rhoj));
  p[i].vxdot -= Aij * wx;
  p[i].vydot -= Aij * wy;
  const double vxi = p[i].vx;
  const double vyi = p[i].vy;
  const double vxj = p[j].vx;
  const double vyj = p[j].vy;
  const double dvxij = vxi - vxj;
  const double dvyij = vyi - vyj;
  p[i].udot += Ci * mj * (dvxij * wx + dvyij * wy);
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
double getRandomJS() {
  return emscripten_random();
}

EMSCRIPTEN_KEEPALIVE
void initializeUniformly(int n,
                         double m,
                         double u, 
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
    particle[i].u = u;
    particle[i].vx = 0.0;
    particle[i].vy = 0.0;
  }
  time = 0.0;
}

EMSCRIPTEN_KEEPALIVE
void zeroVelocity(int n) {
  for (int i = 0; i < n; i++) {
    particle[i].vx = 0.0;
    particle[i].vy = 0.0;
  }
}

EMSCRIPTEN_KEEPALIVE
double getTime(void) {
  return time;
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
double getParticleRho(int i) {
  return particle[i].rho;
}

EMSCRIPTEN_KEEPALIVE
void rebuildTree(int n) {
  double xmin = particle[0].x;
  double xmax = particle[0].x;
  double ymin = particle[0].y;
  double ymax = particle[0].y;

  for (int i = 0; i < n; i++) {
    const double xi = particle[i].x;
    const double yi = particle[i].y;

    qtree.pt[i].w = particle[i].m;
    qtree.pt[i].index = i;
    qtree.pt[i].x = xi;
    qtree.pt[i].y = yi;

    if (xi < xmin) xmin = xi;
    if (xi > xmax) xmax = xi;
    if (yi < ymin) ymin = yi;
    if (yi > ymax) ymax = yi;
  }

  qtree.init(xmin, xmax, ymin, ymax);
  qtree.rebuild(n);
}

EMSCRIPTEN_KEEPALIVE
void setPotentialType(int t) {
  if (t == 0) qtree.setLogPotential();
    else qtree.setInvPotential();
}

EMSCRIPTEN_KEEPALIVE
void computeDensityAndDot(int n, 
                          double G, 
                          double accuracy)
{
  const double gamma_minus_one = gas_gamma - 1.0;

  for (int i = 0; i < n; i++) {
    particle[i].rho = 0.0;
    const tPoint queryi = {particle[i].x, particle[i].y};
    quadtree_box_interact(&qtree.root, i, &queryi, sph_kernel_radius, 
                          &density_summation_callback, reinterpret_cast<void *>(particle));
    particle[i].p = gamma_minus_one * particle[i].u * particle[i].rho;
  }

  for (int i = 0; i < n; i++) {
    particle[i].vxdot = 0.0;
    particle[i].vydot = 0.0;
    particle[i].udot = 0.0;
    const tPoint queryi = {particle[i].x, particle[i].y};
    quadtree_box_interact(&qtree.root, i, &queryi, sph_kernel_radius, 
                          &dot_summation_callback, reinterpret_cast<void *>(particle));
  } 

  if (G <= 0) {
    for (int i = 0; i < n; i++)
      particle[i].phi = 0.0;  
    return;
  }

  const double theta = thetaFromFarness(accuracy);
  for (int i = 0; i < n; i++) {
    const tValueTriad self = qtree.evaluateTarget(particle[i].x, particle[i].y, theta, nearfield_epksq);
    particle[i].phi = G * self.value;
    particle[i].vxdot -= G * self.gradx;
    particle[i].vydot -= G * self.grady;
  }
}

EMSCRIPTEN_KEEPALIVE
void eulerTimestep(int n, double dt) {
  for (int i = 0; i < n; i++) {
    particle[i].u += dt * particle[i].udot;
    particle[i].x += dt * particle[i].vx;
    particle[i].y += dt * particle[i].vy;
    particle[i].vx += dt * particle[i].vxdot;
    particle[i].vy += dt * particle[i].vydot;
  }
  time += dt;
}

EMSCRIPTEN_KEEPALIVE
void boundaryReflection(int n, 
                        double xmin, 
                        double xmax, 
                        double ymin, 
                        double ymax) 
{
  for (int i = 0; i < n; i++) {
    if (particle[i].x < xmin && particle[i].vx < 0.0) particle[i].vx *= -1.0;
    if (particle[i].x > xmax && particle[i].vx > 0.0) particle[i].vx *= -1.0;
    if (particle[i].y < ymin && particle[i].vy < 0.0) particle[i].vy *= -1.0;
    if (particle[i].y > ymax && particle[i].vy > 0.0) particle[i].vy *= -1.0;
  }
}

}
