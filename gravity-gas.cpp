/*
 * WebAssembly simulation of self-gravitating gas in 2D.
 * Uses quadtree for both fixed-radius neighbor search (pressure) and long-range interaction (gravity). 
 * The (ideal) gas is represented with smoothed particle hydrodynamics.
 * The simulation storage elements are: internal energy, gravitational configuration, kinetic energy.
 * There is a viscous loss during compression (as kernels are approaching one another).
 *
 */

#include <emscripten.h>
#include <cmath>
#include <cstring>
#include "qtreehat.hpp"

const double _one_pi = 3.14159265358979323846;

const double gas_gamma = 5.0 / 3.0;

const double nearfield_epksq = 1.0e-6;

const double sph_kernel_h = 5.0;
const double sph_kernel_radius = sph_kernel_h * 2.0;

const double sph_alpha = 1.0;
const double sph_beta = 2.0;
const double sph_eta = 0.01;

// Wendland C2 kernel; support radius = 2 * h
void evaluate_kernel_wc2(double x, 
                         double y, 
                         double h, 
                         double* w,
                         double* wx,
                         double* wy)
{
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
  const double q = std::sqrt(q2);
  const double tmp = 4.0 - 4.0 * q + q2;
  const double W = N * (1.0 + 2.0 * q) * tmp * tmp;
  *w = W;
  const double Wprime = -10.0 * (N / h2) * (2.0 - q) * tmp;
  *wx = Wprime * x;
  *wy = Wprime * y;
  return;
}

const int NMAX = 10000;
const int LEAF = 4;

struct tParticle {
  double m;
  double x;
  double y;
  double u;
  double rho;
  double p;
  double c;
  double vx;
  double vy;
  double phi;
  double udot;
  double vxdot;
  double vydot;
};

struct tDotState {
  double udot;
  double vxdot;
  double vydot;
};

static double time;
static tParticle particle[NMAX];
static tDotState previous_dotstate[NMAX];
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
  tParticle* p = reinterpret_cast<tParticle*>(aux);
  if (i == j) {
    return;
  }
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
  // viscosity
  const double vdotr = dvxij * dxij + dvyij * dyij;
  if (vdotr > 0.0) return;
  const double mu = (sph_kernel_h * vdotr) / (dxij * dxij + dyij * dyij + sph_eta * sph_kernel_h * sph_kernel_h);
  const double cbar = 0.5 * (p[i].c + p[j].c);
  const double rhobar = 0.5 * (rhoi + rhoj);
  const double Bij = mj * (-1.0 * sph_alpha * cbar * mu + sph_beta * mu * mu) / rhobar;
  p[i].vxdot -= Bij * wx;
  p[i].vydot -= Bij * wy;
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
double getRandomJS() {
  return emscripten_random();
}

// Irwin-Hall distribution using n = 12
EMSCRIPTEN_KEEPALIVE
double getApproximateNormal() {
  double sum = 0.0;
  for (int i = 0; i < 12; i++)
    sum += emscripten_random();
  return sum - 6.0;
}

EMSCRIPTEN_KEEPALIVE
double calcLinearMomentumX(int n) {
  double px = 0.0;
  for (int i = 0; i < n; i++)
    px += particle[i].m * particle[i].vx;
  return px;
}

EMSCRIPTEN_KEEPALIVE
double calcLinearMomentumY(int n) {
  double py = 0.0;
  for (int i = 0; i < n; i++)
    py += particle[i].m * particle[i].vy;
  return py;
}

EMSCRIPTEN_KEEPALIVE
double calcAngularMomentumZ(int n) {
  double lz = 0.0;
  for (int i = 0; i < n; i++)
    lz += particle[i].m * (particle[i].x * particle[i].vy - particle[i].y * particle[i].vx);
  return lz;
}

EMSCRIPTEN_KEEPALIVE
double calcTotalEnergy(int n) {
  double K = 0.0;
  double G = 0.0;
  double U = 0.0;
  for (int i = 0; i < n; i++) {
    const double m = particle[i].m;
    K += m * (particle[i].vx * particle[i].vx + particle[i].vy * particle[i].vy);
    G += m * particle[i].phi;
    U += m * particle[i].u; 
  }
  return K / 2.0 + G / 2.0 + U;
}

EMSCRIPTEN_KEEPALIVE
void zeroVelocity(int n) {
  for (int i = 0; i < n; i++) {
    particle[i].vx = 0.0;
    particle[i].vy = 0.0;
  }
}

EMSCRIPTEN_KEEPALIVE
void perturbVelocity(int n, 
                     double sigma) 
{
  for (int i = 0; i < n; i++) {
    particle[i].vx += sigma * getApproximateNormal();
    particle[i].vy += sigma * getApproximateNormal();
  }
}

EMSCRIPTEN_KEEPALIVE
void zeroLinearMomentum(int n) {
  double M = 0.0;
  double px = 0.0;
  double py = 0.0;
  for (int i = 0; i < n; i++) {
    M += particle[i].m;
    px += particle[i].m * particle[i].vx;
    py += particle[i].m * particle[i].vy;
  }
  const double cgvx = px / M;
  const double cgvy = py / M;
  for (int i = 0; i < n; i++) {
    particle[i].vx -= cgvx;
    particle[i].vy -= cgvy;
  }
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
void initializeDisc(int n,
                    double m,
                    double u, 
                    double cx, 
                    double cy, 
                    double R)
{
  int k = 0;
  const double R2 = R * R;
  while (k < n) {
    const double xk = R * (2.0 * getRandomJS() - 1.0);
    const double yk = R * (2.0 * getRandomJS() - 1.0);
    const double rk2 = xk * xk + yk * yk;
    if (rk2 > R2) continue;
    particle[k].x = cx + xk;
    particle[k].y = cy + yk;
    particle[k].m = m;
    particle[k].u = u;
    particle[k].vx = 0.0;
    particle[k].vy = 0.0;
    k++;
  }
  time = 0.0;
}

EMSCRIPTEN_KEEPALIVE
void initializeThreeBody(int n,
                         double m,
                         double u,
                         double cx,
                         double cy,
                         double R,
                         double V)
{
  // Generate 3x cloud centers
  const double theta1 = 2.0 * _one_pi * getRandomJS();
  const double rad1 = getRandomJS() * R;
  const double x1 = cx + std::cos(theta1) * rad1;
  const double y1 = cy + std::sin(theta1) * rad1;

  const double vx1 = getApproximateNormal() * V;
  const double vy1 = getApproximateNormal() * V;

  const double theta2 = 2.0 * _one_pi * getRandomJS();
  const double rad2 = getRandomJS() * R;
  const double x2 = cx + std::cos(theta2) * rad2;
  const double y2 = cy + std::sin(theta2) * rad2;

  const double vx2 = getApproximateNormal() * V;
  const double vy2 = getApproximateNormal() * V;

  const double theta3 = 2.0 * _one_pi * getRandomJS();
  const double rad3 = getRandomJS() * R;
  const double x3 = cx + std::cos(theta3) * rad3;
  const double y3 = cy + std::sin(theta3) * rad3;

  const double vx3 = getApproximateNormal() * V;
  const double vy3 = getApproximateNormal() * V;

  for (int i = 0; i < n; i++) {
    double dx = getApproximateNormal() * R / 8.0;
    double dy = getApproximateNormal() * R / 8.0;

    double dvx = getApproximateNormal() * V / 20.0;
    double dvy = getApproximateNormal() * V / 20.0;

    const double v = getRandomJS();
    if (v < 1.0 / 3.0) {
      dx += x1;
      dy += y1;
      dvx += vx1;
      dvy += vy1;
    } else if (v > 1.0 / 3.0 && v < 2.0 / 3.0) {
      dx += x2;
      dy += y2;
      dvx += vx2;
      dvy += vy2;
    } else {
      dx += x3;
      dy += y3;
      dvx += vx3;
      dvy += vy3;
    }

    particle[i].x = dx;
    particle[i].y = dy;
    particle[i].m = m;
    particle[i].u = u;
    particle[i].vx = dvx;
    particle[i].vy = dvy;
  }

  zeroLinearMomentum(n);
}

EMSCRIPTEN_KEEPALIVE
void forceSpin(int n, 
               double omega,
               bool add)
{
  double sumx = 0.0;
  double sumy = 0.0;
  double summ = 0.0;
  for (int i = 0; i < n; i++) {
    summ += particle[i].m;
    sumx += particle[i].m * particle[i].x;
    sumy += particle[i].m * particle[i].y;
  }
  const double cgx = sumx / summ;
  const double cgy = sumy / summ;
  for (int i = 0; i < n; i++) {
    const double rx = particle[i].x - cgx;
    const double ry = particle[i].y - cgy;
    if (add) {
      particle[i].vx -= omega * ry;
      particle[i].vy += omega * rx;
    } else {
      particle[i].vx = -1.0 * omega * ry;
      particle[i].vy = omega * rx;
    }
  }
}

EMSCRIPTEN_KEEPALIVE
void multiplySpecificEnergy(int n,
                            double s)
{
  for (int i = 0; i < n; i++) {
    particle[i].u *= s;
  }
}

EMSCRIPTEN_KEEPALIVE
double getTime(void) {
  return time;
}

EMSCRIPTEN_KEEPALIVE
void updateTime(double dt) {
  time += dt;
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
double averageLeafsize(void) {
  return qtree.average_leafsize();
}

EMSCRIPTEN_KEEPALIVE
double averageDepth(void) {
  return qtree.average_depth();
}

EMSCRIPTEN_KEEPALIVE
bool addParticleAt(int n, 
                   double m, 
                   double u,
                   double x,
                   double y,
                   double vx, 
                   double vy) 
{
  if (n >= NMAX) return false;
  particle[n].x = x;
  particle[n].y = y;
  particle[n].vx = vx;
  particle[n].vy = vy;
  particle[n].u = u;
  particle[n].m = m;
  return true;
}

EMSCRIPTEN_KEEPALIVE
bool addGaussianAt(int n,
                   int g, 
                   double m, 
                   double u,
                   double x,
                   double y,
                   double sigma,
                   double vx, 
                   double vy) 
{
  if (n + g >= NMAX) return false;
  for (int i = n; i < n + g; i++) {
    particle[i].x = x + sigma * getApproximateNormal();
    particle[i].y = y + sigma * getApproximateNormal();
    particle[i].vx = vx;
    particle[i].vy = vy;
    particle[i].u = u;
    particle[i].m = m;
  }
  return true;
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
  const double sqrt_g_gm1 = std::sqrt(gas_gamma * gamma_minus_one);

  for (int i = 0; i < n; i++) {
    particle[i].rho = 0.0;
    const tPoint queryi = {particle[i].x, particle[i].y};
    quadtree_box_interact(&qtree.root, i, &queryi, sph_kernel_radius, 
                          &density_summation_callback, reinterpret_cast<void *>(particle));
    particle[i].p = gamma_minus_one * particle[i].u * particle[i].rho;
    particle[i].c = sqrt_g_gm1 * std::sqrt(particle[i].u);
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
void eulerUpdate(int n, double dt) {
  for (int i = 0; i < n; i++) {
    particle[i].u += dt * particle[i].udot;
    particle[i].x += dt * particle[i].vx;
    particle[i].y += dt * particle[i].vy;
    particle[i].vx += dt * particle[i].vxdot;
    particle[i].vy += dt * particle[i].vydot;
  }
}

EMSCRIPTEN_KEEPALIVE
void saveDotState(int n) {
  for (int i = 0; i < n; i++) {
    previous_dotstate[i].udot = particle[i].udot;
    previous_dotstate[i].vxdot = particle[i].vxdot;
    previous_dotstate[i].vydot = particle[i].vydot;
  }
}

EMSCRIPTEN_KEEPALIVE
void firstUpdate(int n, double dt) {
  const double A = dt * dt / 2.0;
  for (int i = 0; i < n; i++) {
    particle[i].x += dt * particle[i].vx + A * particle[i].vxdot;
    particle[i].y += dt * particle[i].vy + A * particle[i].vydot;
    //particle[i].u += dt * particle[i].udot;
  }
}

EMSCRIPTEN_KEEPALIVE
void secondUpdate(int n, double dt) {
  const double A = dt / 2.0;
  for (int i = 0; i < n; i++) {
    particle[i].vx += A * (particle[i].vxdot + previous_dotstate[i].vxdot);
    particle[i].vy += A * (particle[i].vydot + previous_dotstate[i].vydot);
    //particle[i].u  += A * (-1.0 * previous_dotstate[i].udot + particle[i].udot);
    particle[i].u  += A * (particle[i].udot + previous_dotstate[i].udot);
  }
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
