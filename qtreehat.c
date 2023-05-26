/*
 * qtreehat.c
 *
 * All-to-all interaction evaluation using a point-region
 * quad-tree construction. Matrix-free and mesh-free.
 * Optionally provide target points (t) separate from sources (p).
 *
 * USAGE: wp = qtreehat(xp, yp, vp, maxleaf, ep, epk, pot);
 *        wt = qtreehat(xp, yp, vp, maxleaf, ep, epk, pot, xt, yt);
 *
 *   wp(i) = sum_j { vp(j) * log(|r(i) - r(j)|) },    if pot = 'log'
 *   wp(i) = sum_j {-vp(j) / |r(i) - r(j)| },         if pot = 'inv'
 *
 * Set ep <= 0 to compute the exact reference sum (brute force).
 *
 * Set epk > 0 to (roughly) regularize the (very) near field evaluation.
 * If epk = 0, the self-terms are dropped.
 *
 * (for a gaussian kernel with bandwidth h, epk = 0.7493 * h, for log-potential)
 *
 * Also calculates the gradient of wp(i) w.r.t. r(i)_x, and r(i)_y.
 * So wp has three columns: {value, grad_x, grad_y}
 *
 * Two outputs: [wp, nno] = qptreehat(..), then nno = number of nodes in the tree
 * Zero outputs: wp will be stored in "ans" + debug/stats print-outs are shown
 *
 * COMPILE: build_mex('qtreehat');
 * TEST:    test_qtreehat.m
 *
 */

#include "mex.h"
#include <stdint.h>
#include <math.h>
#include <string.h>

#define THEMALLOC mxMalloc
#define THEFREE   mxFree
#define THEPRINTF mexPrintf

#include "qtreehat.h"

#define THEERRMSG mexErrMsgTxt
#define THEWRNMSG mexWarnMsgTxt

/* ------------------------------------------------------------ */
/* Some MEX I/O helpers                                         */
/* ------------------------------------------------------------ */

bool isDoubleRealMatrix(const mxArray* a) {
  if (a == NULL) return false;
  if (mxIsEmpty(a)) return false;
  if (!mxIsDouble(a)) return false;
  if (mxGetNumberOfDimensions(a) != 2) return false;
  if (mxIsComplex(a)) return false;
  return true;
}

bool isDoubleRealVector(const mxArray* a) {
  if (!isDoubleRealMatrix(a)) return false;
  int ma = mxGetM(a);
  int na = mxGetN(a);
  if (ma != 1 && na != 1) return false;
  return true;
}

/* ------------------------------------------------------------ */
/* MEX main                                                     */
/* ------------------------------------------------------------ */

#define NUM_INARGS   7
#define ARG_XP       prhs[0]
#define ARG_YP       prhs[1]
#define ARG_VP       prhs[2]
#define ARG_MAXLEAF  prhs[3]
#define ARG_EP       prhs[4]
#define ARG_EPK      prhs[5]
#define ARG_POTSTR   prhs[6]

#define NUM_INARGS_T 9
#define ARG_XT       prhs[7]
#define ARG_YT       prhs[8]

#define MAX_OUTARGS  2
#define OUT_WP       plhs[0]
#define OUT_NNODES   plhs[1]

void mexFunction(int nlhs, 
                 mxArray** plhs,
                 int nrhs, 
                 const mxArray** prhs)
{
  if (nrhs != NUM_INARGS && nrhs != NUM_INARGS_T) {
    THEPRINTF("USAGE: wp = qtreehat(xp, yp, vp, maxleaf, ep, epk, pot [,xt, yt]);\n");
    THEERRMSG("Incorrect number of input arguments provided");
  }

  if (nlhs > MAX_OUTARGS) {
    THEERRMSG("Too many output arguments requested.");
  }

  const bool has_separate_targets = (nrhs == NUM_INARGS_T);

  if (!isDoubleRealVector(ARG_XP)) {
    THEERRMSG("xp must be a real vector"); 
  }

  const int nxp = mxGetNumberOfElements(ARG_XP);
  const double* pXP = mxGetPr(ARG_XP);

  if (!isDoubleRealVector(ARG_YP)) {
    THEERRMSG("yp must be a real vector"); 
  }

  const int nyp = mxGetNumberOfElements(ARG_YP);
  const double* pYP = mxGetPr(ARG_YP);

  if (nyp != nxp) {
    THEERRMSG("xp and yp must have the same number of elements");
  }

  if (!isDoubleRealVector(ARG_VP)) {
    THEERRMSG("vp must be a real vector"); 
  }

  const int nvp = mxGetNumberOfElements(ARG_VP);
  const double* pVP = mxGetPr(ARG_VP);

  if (nvp != nxp) {
    THEERRMSG("vp must have the same number of elements as xp, yp"); 
  }

  int nxt = 0;
  const double* pXT = pXP;
  const double* pYT = pYP;
  if (has_separate_targets) {
    if (!isDoubleRealVector(ARG_XT)) {
      THEERRMSG("xt must be a real vector"); 
    }

    nxt = mxGetNumberOfElements(ARG_XT);

    if (!isDoubleRealVector(ARG_YT)) {
      THEERRMSG("yt must be a real vector"); 
    }

    if (nxt != mxGetNumberOfElements(ARG_YT)) {
      THEERRMSG("xt and yt must have same number of elements");
    }

    pXT = mxGetPr(ARG_XT);
    pYT = mxGetPr(ARG_YT);
  }

  if (!(isDoubleRealVector(ARG_EP) && mxGetNumberOfElements(ARG_EP) == 1)) {
    THEERRMSG("ep must be a real-valued scalar"); 
  }

  const double EP = *mxGetPr(ARG_EP);
  const bool EP_is_valid = (EP > 0.0 & EP < 1.0);

  if (!EP_is_valid && nlhs == 2) {
    THEERRMSG("Cannot return a second output when ep is out of range (no quadtree is built)");
  }

  if (!(isDoubleRealVector(ARG_MAXLEAF) && mxGetNumberOfElements(ARG_MAXLEAF) == 1)) {
    THEERRMSG("maxleaf must be a real-valued scalar"); 
  }

  const int maxLeaf = (int) floor(*mxGetPr(ARG_MAXLEAF));

  if (maxLeaf <= 0) {
    THEERRMSG("maxLeaf >= 1 required");
  }

  if (!(isDoubleRealVector(ARG_EPK) && mxGetNumberOfElements(ARG_EPK) == 1)) {
    THEERRMSG("epk must be a real-valued scalar"); 
  }

  const double EPK = *mxGetPr(ARG_EPK);
  if (EPK < 0.0) {
    THEERRMSG("epk >= 0 required");
  }

  if (!mxIsChar(ARG_POTSTR)) {
    THEERRMSG("pot must be a string");
  }

  const char* pot_str = mxArrayToString(ARG_POTSTR);
  const bool use_log_pot = (strcmp(pot_str, "log") == 0);
  const bool use_inv_pot = (strcmp(pot_str, "inv") == 0);

  if (!use_inv_pot && !use_log_pot) {
    THEERRMSG("pot must be either \"log\" or \"inv\" (case sensitive)");
  }

  const double epksq = EPK * EPK;

  PotentialFuncPtr pfunc  = (use_log_pot ? &logr_potential_ : &invr_potential_);
  QuadrupoleFuncPtr qfunc = (use_log_pot ? &eval_logr_quadrupole_ : &eval_invr_quadrupole_);

  const int Np = nxp;
  const int Nt = (has_separate_targets ? nxt : nxp);
  OUT_WP = mxCreateDoubleMatrix(Nt, 3, mxREAL);
  double* pWP = mxGetPr(OUT_WP);

  if (!EP_is_valid) {
    for (int i = 0; i < Nt; i++) {
      const double xi = pXT[i];
      const double yi = pYT[i];
      double sumf = 0.0;
      double sumfx = 0.0;
      double sumfy = 0.0;
      for (int j = 0; j < Np; j++) {
        const double xj = pXP[j];
        const double yj = pYP[j];
        const tValueTriad vt = (*pfunc)(xi, yi, xj, yj, epksq);  // gradient w.r.t. 1st tuple (i)
        sumf  += vt.value * pVP[j];
        sumfx += vt.gradx * pVP[j];
        sumfy += vt.grady * pVP[j];
      }
      pWP[i] = sumf;
      pWP[i + Nt] = sumfx;
      pWP[i + 2 * Nt] = sumfy;
    }

    return;
  }

  const int maxNumberOfNodes = estimateMaxNumNodes(Np, maxLeaf);

  tPointPayload* pt = THEMALLOC(sizeof(tPointPayload) * Np);
  tPointPayload* pt_scratch = THEMALLOC(sizeof(tPointPayload) * Np);
  tQuadTree* nodes = THEMALLOC(sizeof(tQuadTree) * maxNumberOfNodes);
  
  double xmin = pXP[0];
  double xmax = pXP[0];
  double ymin = pYP[0];
  double ymax = pYP[0];
  for (int i = 0; i < Np; i++) {
    pt[i].x = pXP[i];
    pt[i].y = pYP[i];
    pt[i].w = pVP[i]; // load up weight value for each point also
    xmin = (pXP[i] < xmin ? pXP[i] : xmin);
    xmax = (pXP[i] > xmax ? pXP[i] : xmax);
    ymin = (pYP[i] < ymin ? pYP[i] : ymin);
    ymax = (pYP[i] > ymax ? pYP[i] : ymax);
  }

  const double hbwx = (xmax - xmin) / 2.0;
  const double hbwy = (ymax - ymin) / 2.0;

  tQuadTree rootNode;
  rootNode.cb.x = (xmin + xmax) / 2.0;
  rootNode.cb.y = (ymin + ymax) / 2.0;
  rootNode.hbw = (hbwx > hbwy ? hbwx : hbwy);
  rootNode.hbw *= (1.0 + 1.0e-10);

  const int maxDepth = 50;
  const int maxInLeaf = maxLeaf;
  const int rootLevel = 0;

  int nno = build_quadtree(&rootNode, 
                           maxInLeaf, 
                           maxDepth, 
                           rootLevel, 
                           Np, 
                           pt, 
                           pt_scratch, 
                           maxNumberOfNodes, 
                           nodes);

  if (nlhs == 0) {
    THEPRINTF("[%s]: maxLeaf = %i\n", __func__, maxLeaf);
    THEPRINTF("[%s]: nno = %i (allocated = %i)\n", __func__, nno, maxNumberOfNodes);

    const int nn_in_tree = count_quadtree_nodes(&rootNode);
    THEPRINTF("[%s]: nn_in_tree = %i\n", __func__, nn_in_tree);

    const int n_in_tree = count_quadtree_points(&rootNode);
    THEPRINTF("[%s]: n_in_tree = %i\n", __func__, n_in_tree);

    int total_box_count = quadtree_box_query_count(&rootNode, &(rootNode.cb), rootNode.hbw);
    THEPRINTF("[%s]: n_in_tree (total box query) = %i\n", __func__, total_box_count);

    const int max_depth = count_maximum_depth(&rootNode, 0);
    THEPRINTF("[%s]: max_depth = %i (zero is root)\n", __func__, max_depth);

    const double avg_depth = count_average_depth(&rootNode, 0.0);
    THEPRINTF("[%s]: avg_depth = %f (zero is root)\n", __func__, avg_depth);

    tMomentStats refstats;
    calc_source_summary(&refstats, Np, pt, rootNode.cb.x, rootNode.cb.y);
    compare_summaries(&refstats, &(rootNode.stats));
  }

  const double theta = EP * EP;
  for (int i = 0; i < Nt; i++) {
    const double xi = pXT[i];
    const double yi = pYT[i];
    const tPoint ith_query = {xi, yi};
    const tValueTriad wi = quadtree_sum_at_point(&rootNode, 
                                                 &ith_query, 
                                                 theta, 
                                                 epksq, 
                                                 pfunc, 
                                                 qfunc);
    pWP[i] = wi.value;
    pWP[i + Nt] = wi.gradx;
    pWP[i + 2 * Nt] = wi.grady;
  }

  THEFREE(pt);
  THEFREE(pt_scratch);
  THEFREE(nodes);

  if (nlhs == 2) {
    OUT_NNODES = mxCreateDoubleScalar((double)(nno + 1));
  }

  return;
}
