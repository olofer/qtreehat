// BUILD: g++ -std=c++14 -Wall -o test-index++.exe -O2 test-index.cpp

// TODO: verify that the tree works for variable number of points N, as long as N <= NMAX @ fixed allocation

#include <iostream>
#include <random>

#include <cstring>  // for std::memset()

#define _QTREEHAT_NOINDEX_
#define _QTREEHAT_NOPRINTF_
#include "qtreehat.h"

constexpr int estMaxNodes(int N, int L)
{
  if (L <= 0) L = 1;
  const double f1 = 16.0 / 3.0;
  const double f2 = (double) (L + 1);
  const double f = (f1 < f2 ? f1 : f2);
  const double n_ovr_l = ((double) N) / L;
  return (int) (f * n_ovr_l) + 128;
}

template <int N, int L, int M = estMaxNodes(N, L)>
struct tQuadTreeHandle
{
  tPointPayload pt[N];
  tPointPayload pt_scratch[N];
  tQuadTree nodes_store[M];
  tQuadTree root;

  PotentialFuncPtr pfunc;
  QuadrupoleFuncPtr qfunc;

  void setLogPotential() {
    pfunc = &logr_potential_;
    qfunc = &eval_logr_quadrupole_;
  }

  void setInvPotential() {
    pfunc = &invr_potential_;
    qfunc = &eval_invr_quadrupole_;
  }

  int capacityNodes() const { return M; }
  int capacityPoints() const { return N; }
  int capacityLeaf() const { return L; }

  bool init(double xmin, 
            double xmax, 
            double ymin, 
            double ymax) 
  {
    if (xmin >= xmax || ymin >= ymax) return false;

    const double hbwx = (xmax - xmin) / 2.0;
    const double hbwy = (ymax - ymin) / 2.0;
    const double cbx = (xmin + xmax) / 2.0;
    const double cby = (ymin + ymax) / 2.0;
    const double eps_mult = 1.0e-10;

    root.cb.x = cbx;
    root.cb.y = cby;
    root.hbw = (hbwx > hbwy ? hbwx : hbwy);
    root.hbw *= (1.0 + eps_mult);

    return true;
  }

  int rebuild(int numpts, 
              int maxDepth = 50) 
  {
    if (numpts <= 0 || numpts > N) return 0;

    const int rootLevel = 0;
    const int nno = build_quadtree(&root, 
                                   L, 
                                   maxDepth, 
                                   rootLevel, 
                                   numpts, 
                                   pt, 
                                   pt_scratch, 
                                   M, 
                                   nodes_store);
    return nno;
  }

  int maximum_leafsize() const {
    return count_maximum_in_leaf(&root);
  }

  double average_leafsize() const {
    return count_average_in_leaf(&root);
  }

  int maximum_depth() const {
    return count_maximum_depth(&root, 0);
  }

  double average_depth() const {
    return count_average_depth(&root, 0);
  }

  int countPoints() const {
    return count_quadtree_points(&root);
  }

  int countNodes() const {
    return count_quadtree_nodes(&root);
  }

  tValueTriad evaluateTarget(double x, 
                             double y,
                             double theta,
                             double epksq) const 
  {
    const tPoint ith_query = {x, y};
    const tValueTriad triad = quadtree_sum_at_point(&root, 
                                                    &ith_query, 
                                                    theta, 
                                                    epksq, 
                                                    pfunc, 
                                                    qfunc);
    return triad;
  }
};

int main(int argc, 
         const char** argv)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> Normal(0.0, 1.0);
  std::uniform_real_distribution<double> Uniform(0.0, 1.0);

  if (argc != 4) {
    std::cout << "usage: " << argv[0] << " numpts D qhw" << std::endl;
    return 0;
  }

  const int numpts = std::atoi(argv[1]);
  const double D = std::atof(argv[2]);    // domain width
  const double qhw = std::atof(argv[3]);  // query half-width

  const int NMAX = 25000;
  const int LEAF = 8;

  if (numpts <= 0 || numpts > NMAX || D <= 0.0 || qhw <= 0.0) {
    std::cout << "infeasible argments" << std::endl;
    return 1;
  }

  tQuadTreeHandle<NMAX, LEAF> qtHandle;

  for (int i = 0; i < numpts; i++) {
    qtHandle.pt[i].index = i;
    qtHandle.pt[i].x = D * Uniform(gen) - D / 2.0;
    qtHandle.pt[i].y = D * Uniform(gen) - D / 2.0;
    qtHandle.pt[i].w = Normal(gen);

    if (i == numpts - 1) {
      std::cout << "w[" << i << "] = " << qtHandle.pt[i].w << std::endl;
    }
  }

  const double xmin = -D / 2.0;
  const double xmax = +D / 2.0;
  const double ymin = -D / 2.0;
  const double ymax = +D / 2.0;

  qtHandle.init(xmin, xmax, ymin, ymax);
  qtHandle.rebuild(numpts);

  std::cout << "    points = " << qtHandle.countPoints() << " (capacity = " << qtHandle.capacityPoints() << ")" << std::endl;
  std::cout << "    nodes  = " << qtHandle.countNodes()  << " (capacity = " << qtHandle.capacityNodes() << ")" << std::endl;
  std::cout << "<leafsize> = " << qtHandle.average_leafsize() << " (maximum = " << qtHandle.maximum_leafsize() << ")" << std::endl;
  std::cout << "   <depth> = " << qtHandle.average_depth()    << " (maximum = " << qtHandle.maximum_depth() << ")" << std::endl;

  // TODO: test a subsample of target evaluations ...
  qtHandle.setLogPotential();

  // TODO: run the test query for correctness ...

  return 0;
}
