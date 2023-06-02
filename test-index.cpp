/*
 * Basic test of C++ interface for fixed-maximum-size (by template parameters) treecode.
 * Similar to test-index.c but goes a little further and also tests the target point eval() member.
 * 
 * BUILD: g++ -std=c++14 -Wall -o test-index++.exe -O2 test-index.cpp
 *        clang++ -std=c++14 -Wall -o test-index++.exe -O2 test-index.cpp
 */


#include <iostream>
#include <random>

#include <cmath>
#include <cstring>  // for std::memset()
#include <vector>

#include "qtreehat.hpp"

void enumerateFunction(int i, 
                       int j, 
                       void* aux)
{
  if (aux == nullptr) return;
  uint32_t* histo = (uint32_t *) aux;
  histo[j]++;
}

void count_all_inside(const tPointPayload* pt, 
                      int n,
                      std::vector<uint32_t>& histo,
                      double xmin,
                      double xmax,
                      double ymin,
                      double ymax)
{
  for (int i = 0; i < n; i++) {
    const double xi = pt[i].x;
    const double yi = pt[i].y;
    const bool x_inside = (xi >= xmin && xi < xmax);
    const bool y_inside = (yi >= ymin && yi < ymax);
    if (x_inside && y_inside) {
      histo[i]++;
    }
  }
}

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

  std::cout << "template parameters N, L = " << NMAX << ", " << LEAF << std::endl;

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

  // Run a brute force test of nearest-neighbor detection  
  std::vector<uint32_t> H0(numpts, 0);
  std::vector<uint32_t> H1(numpts, 0);

  std::cout << "brute force testing of interaction search box ..." << std::endl;

  qtHandle.restore_scratch_in_original_order(numpts);

  const double query_box_halfwidth = qhw;

  for (int i = 0; i < numpts; i++) {
    const tPoint query_pt_i = {qtHandle.pt_scratch[i].x, qtHandle.pt_scratch[i].y};
    const int index_query_i = i;

    quadtree_box_interact(&qtHandle.root,
                          index_query_i, 
                          &query_pt_i, 
                          query_box_halfwidth,
                          &enumerateFunction,
                          reinterpret_cast<void *>(H0.data()));

    count_all_inside(qtHandle.pt_scratch, 
                     numpts, 
                     H1, 
                     query_pt_i.x - query_box_halfwidth, 
                     query_pt_i.x + query_box_halfwidth,
                     query_pt_i.y - query_box_halfwidth, 
                     query_pt_i.y + query_box_halfwidth);

    const bool is_still_equal = (memcmp(H0.data(), H1.data(), sizeof(uint32_t) * numpts) == 0);

    if (!is_still_equal) {
      std::cout << "interaction error @ " << i << std::endl;
      return -1;
    }
  }

  uint32_t hmin = 1000000000;
  uint32_t hmax = 0;
  for (int i = 0; i < numpts; i++) {
    if (H0[i] < 1) {
      std::cout << "never self interacted @ " << i << std::endl;
      return -2;
    }
    if (H0[i] < hmin) hmin = H0[i];
    if (H0[i] > hmax) hmax = H0[i];
  }
  std::cout << "histogram min, max = " << hmin << ", " << hmax << std::endl;

  // Check accuracy of tree code evaluation as controlled by the theta parameter

  std::cout << "brute force testing target point evaluation accuracies ..." << std::endl;
  qtHandle.setLogPotential();

  std::vector<tValueTriad> referenceValues(numpts);

  const double theta_ref = 1.0e-12;
  const double epksq = 1.0e-8;
  double reference_sum = 0.0; 

  for (int i = 0; i < numpts; i++)
  {
    referenceValues[i] = qtHandle.evaluateTarget(qtHandle.pt_scratch[i].x, 
                                                 qtHandle.pt_scratch[i].y, 
                                                 theta_ref, 
                                                 epksq);
    
    reference_sum += std::fabs(referenceValues[i].value);
  }

  std::vector<tValueTriad> testValues(numpts);
  const std::vector<double> ep_tests = {1.0 / 256, 1.0 / 128, 1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8, 1.0 / 4, 1.0 / 2};

  // Recall that "ep" can be regarded as the (region size) / (distance to region) ratio; and "theta" is its square

  for (double ep : ep_tests)
  {
    const double theta = thetaFromFarness(ep); // ep * ep;
    std::cout << "theta = " << theta << " (ep = " << ep << ")";

    double error_sum = 0.0;
    for (int i = 0; i < numpts; i++)
    {
      testValues[i] = qtHandle.evaluateTarget(qtHandle.pt_scratch[i].x, 
                                              qtHandle.pt_scratch[i].y, 
                                              theta, 
                                              epksq);
      
      error_sum += std::fabs(testValues[i].value - referenceValues[i].value);
    }

    std::cout << "; error = " << error_sum / reference_sum << std::endl;
  }

  return 0;
}
