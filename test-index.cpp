// BUILD: g++ -std=c++14 -Wall -o test-index++.exe -O2 test-index.cpp

// TODO: verify that the tree works for variable number of points N, as long as N <= NMAX @ fixed allocation

#include <iostream>
#include <random>

#include <cstring>  // for std::memset()

#define _QTREEHAT_NOINDEX_
#define _QTREEHAT_NOPRINTF_
#include "qtreehat.h"

// TODO: make a C++ template handle for the quadtree
// probably also <cmath> will be required, "constexpr"

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

  bool init(double xmin, 
            double xmax, 
            double ymin, 
            double ymax) 
  {
    // TODO: add init box method etc.. 
    return false;
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

  // TODO: and then verify that the standard C callbacks can be applied 
};

int main(int argc, 
         const char** argv)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> Normal(0.0, 1.0);
  std::uniform_real_distribution<double> Uniform(0.0, 1.0);

  std::cout << "usage: " << argv[0] << " (blah)" << std::endl;
  // TODO: parse inputs n D d ?

  const int NMAX = 25000;
  const int LEAF = 8;

  tQuadTreeHandle<NMAX, LEAF> qtHandle;

  for (int i = 0; i < NMAX; i++) {
    qtHandle.pt[i].index = i;
    qtHandle.pt[i].x = Uniform(gen);
    qtHandle.pt[i].y = Uniform(gen);
    qtHandle.pt[i].w = Normal(gen);

    if (i == NMAX - 1) {
      std::cout << qtHandle.pt[i].w << std::endl;
    }
  }

  qtHandle.rebuild(0);

  return 0;
}
