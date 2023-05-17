/*
 * qtreehat.c
 *
 * All-to-all interaction evaluation using a point-region
 * quad-tree construction. Matrix-free and mesh-free.
 *
 * USAGE: wp = qtreehat(xp, yp, vp, ep, maxleaf, epk, pot);
 *
 *   wp(i) = sum_j { vp(j) * log(|r(i) - r(j)|) },    if pot = 'log'
 *   wp(i) = sum_j {-vp(j) / |r(i) - r(j)| },         if pot = 'inv'
 *
 * Set ep <= 0 to compute the exact reference sum.
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

/*
  Each node stores sums of:

    0) v
    1) v * [x; y]
    2) v * [x; y] * [x; y]'

  This is information required to compute quadrupole approximations from the sources under a node.
  The sums are built up hierarchically; children nodes are merged using translation rules + direct sums.
*/

#include "mex.h"
#include <stdint.h>
#include <math.h>
#include <string.h>

#define THEMALLOC mxMalloc
#define THEFREE   mxFree
#define THEPRINTF mexPrintf
#define THEERRMSG mexErrMsgTxt
#define THEWRNMSG mexWarnMsgTxt

typedef struct tValueTriad {
  double value;
  double gradx;
  double grady;
} tValueTriad;

// This is directly computed for leaf nodes;
// non-leaf nodes compute this by (correctly) merging the stats from its 4 child nodes.
typedef struct tMomentStats {
  double w;
  double wx;   // "dipole" terms
  double wy;
  double wxx;  // "quadrupole" terms
  double wxy;
  double wyy;
  int nw;
} tMomentStats;

typedef tValueTriad (*PotentialFuncPtr)(double, double, double, double, double);
typedef tValueTriad (*QuadrupoleFuncPtr)(double, double, const tMomentStats*);

// potential: log(|R|)
// gradient: w.r.t. (x1, y1)
tValueTriad logr_potential_(double x1, 
                            double y1, 
                            double x2, 
                            double y2,
                            double epsq)
{
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  const double rsq = dx * dx + dy * dy + epsq;
  tValueTriad vt = {0.0, 0.0, 0.0};
  if (rsq == 0.0) return vt;
  vt.value = log(rsq) / 2.0;
  vt.gradx = dx / rsq;
  vt.grady = dy / rsq;
  return vt;
}

// potential: -1/|R|
// gradient: w.r.t. (x1, y1)
tValueTriad invr_potential_(double x1, 
                            double y1, 
                            double x2, 
                            double y2,
                            double epsq)
{
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  const double rsq = dx * dx + dy * dy + epsq;
  tValueTriad vt = {0.0, 0.0, 0.0};
  if (rsq == 0.0) return vt;
  vt.value = -1.0 / sqrt(rsq);
  vt.gradx = -dx * vt.value / rsq;  // same as: dx / rsq^(3/2)
  vt.grady = -dy * vt.value / rsq;  //          dy / rsq^(3/2)
  return vt;
}

// Pre-specify max-nodes and max-leaves; pre-allocate sufficient space?
// bounding box = [xmin, xmax, ymin, ymax] = [cx-hw, cx+hw, cy-hw, cy+hw]

typedef struct tPoint {
  double x;
  double y;
} tPoint;

typedef struct tPointPayload {
  double x;
  double y;
  double w;
} tPointPayload;

typedef struct tQuadTree {
  tPoint cb;   // box center = (cx, cy)
  double hbw;  // half box width, box is a square always
  tMomentStats stats;

  int npts;    // if leaf node: num. points in leaf, otherwise sum-total of points in children
  tPointPayload* p;   // NULL for non-leaf nodes; otherwise actual point data storage

  struct tQuadTree* pp; // child-nodes, one per quadrant
  struct tQuadTree* pm; // NULL if empty
  struct tQuadTree* mm;
  struct tQuadTree* mp;
} tQuadTree;

tValueTriad eval_logr_quadrupole_(double X, 
                                  double Y, 
                                  const tMomentStats* m)
{
  const double Q = m->w;
  const double Qx = m->wx;
  const double Qy = m->wy;
  const double Qxx = m->wxx;
  const double Qxy = m->wxy;
  const double Qyy = m->wyy;

  const double R2 = X * X + Y * Y;
  const double f0 = Q * log(R2) / 2.0;
  const double f1 = -1.0 * (X * Qx + Y * Qy) / R2;
  const double Cxx = - X * X + Y * Y;
  const double Cyy = - Cxx; // X * X - Y * Y;
  const double Cxy = -2.0 * X * Y;
  const double R4 = R2 * R2;
  const double f2 = 0.5 * (Cxx * Qxx + Cxy * Qxy + Cyy * Qyy) / R4;

  const double g0x = Q * X / R2;
  const double g0y = Q * Y / R2;
  const double g1x = - (Cxx * Qx + Cxy * Qy) / R4;
  const double g1y = - (Cxy * Qx + Cyy * Qy) / R4;

  const double R6 = R4 * R2;
  const double CX = 2.0 * X * (X * X - 3.0 * Y * Y);
  const double CY = 2.0 * Y * (Y * Y - 3.0 * X * X);
  const double g2x = 0.5 * (CX * (Qxx - Qyy) - CY * Qxy) / R6;
  const double g2y = 0.5 * (CY * (Qyy - Qxx) - CX * Qxy) / R6;

  const tValueTriad vt = {f0 + f1 + f2, g0x + g1x + g2x, g0y + g1y + g2y};
  return vt;
}

tValueTriad eval_invr_quadrupole_(double X, 
                                  double Y, 
                                  const tMomentStats* m)
{
  const double Q = m->w;
  const double Qx = m->wx;
  const double Qy = m->wy;
  const double Qxx = m->wxx;
  const double Qxy = m->wxy;
  const double Qyy = m->wyy;

  const double R2 = X * X + Y * Y;
  const double R = sqrt(R2);
  const double Xhat = X / R;
  const double Yhat = Y / R;

  const double f0 = -1.0 * Q / R;
  const double f1 = -1.0 * (Qx * Xhat + Qy * Yhat) / R2;

  const double XXh = Xhat * Xhat;
  const double XYh = Xhat * Yhat;
  const double YYh = Yhat * Yhat;

  const double R3 = R * R2;
  const double f2 = (Qxx * (0.5 * YYh - XXh) + Qxy * (-3.0 * XYh) + Qyy * (0.5 * XXh - YYh)) / R3;

  const double g0x = Q * Xhat / R2;
  const double g0y = Q * Yhat / R2;

  const double CX_x = 2.0 * XXh - YYh;
  const double CX_y = 3.0 * XYh; // same as CY_x
  const double CY_y = 2.0 * YYh - XXh;

  const double g1x = (CX_x * Qx + CX_y * Qy) / R3;
  const double g1y = (CX_y * Qx + CY_y * Qy) / R3;

  const double CX_xx = 3.0 * Xhat * (XXh - 1.5 * YYh);
  const double CY_xx = 3.0 * Yhat * (2.0 * XXh - 0.5 * YYh);

  const double CX_xy = -3.0 * Yhat * (YYh - 4.0 * XXh);
  const double CY_xy = -3.0 * Xhat * (XXh - 4.0 * YYh);

  const double CX_yy = 3.0 * Xhat * (2.0 * YYh - 0.5 * XXh);
  const double CY_yy = 3.0 * Yhat * (YYh - 1.5 * XXh);

  const double R4 = R2 * R2;
  const double g2x = (Qxx * CX_xx + Qxy * CX_xy + Qyy * CX_yy) / R4;
  const double g2y = (Qxx * CY_xx + Qxy * CY_xy + Qyy * CY_yy) / R4;

  const tValueTriad vt = {f0 + f1 + f2, g0x + g1x + g2x, g0y + g1y + g2y};
  return vt;
}

// Up-to-2nd order moments of the set of source points w.r.t to (refx, refy).
void calc_source_summary(tMomentStats* m,
                         int npts, 
                         const tPointPayload* p,
                         double refx,
                         double refy) 
{
  memset(m, 0, sizeof(tMomentStats));
  for (int i = 0; i < npts; i++) {
    const double xi = p[i].x - refx;
    const double yi = p[i].y - refy;
    const double wi = p[i].w;
    m->w   += wi;
    m->wx  += wi * xi;
    m->wy  += wi * yi;
    m->wxx += wi * xi * xi;
    m->wxy += wi * xi * yi;
    m->wyy += wi * yi * yi;
  }
  m->nw += npts;
  return;
}

// Use this to aggregate (translated) child node stats
void add_source_translation(tMomentStats* m,
                            const tMomentStats* msrc,
                            double bx,
                            double by)
{
  m->w   += (msrc->w);
  m->wx  += (msrc->wx) + (msrc->w) * bx;
  m->wy  += (msrc->wy) + (msrc->w) * by;
  m->wxx += (msrc->wxx) + 2.0 * bx * (msrc->wx) + (msrc->w) * bx * bx;
  m->wxy += (msrc->wxy) + bx * (msrc->wy) + by * (msrc->wx) + (msrc->w) * bx * by;
  m->wyy += (msrc->wyy) + 2.0 * by * (msrc->wy) + (msrc->w) * by * by;
  m->nw  += (msrc->nw);
  return;
}

void compare_summaries(const tMomentStats* S1,
                       const tMomentStats* S2)
{
  THEPRINTF("[%s]: S1.nw,  S2.nw  = %i, %i\n", __func__, S1->nw, S2->nw);
  THEPRINTF("[%s]: S1.w,   S2.w   = %e, %e\n", __func__, S1->w, S2->w);
  THEPRINTF("[%s]: S1.wx,  S2.wx  = %e, %e\n", __func__, S1->wx, S2->wx);
  THEPRINTF("[%s]: S1.wy,  S2.wy  = %e, %e\n", __func__, S1->wy, S2->wy);
  THEPRINTF("[%s]: S1.wxx, S2.wxx = %e, %e\n", __func__, S1->wxx, S2->wxx);
  THEPRINTF("[%s]: S1.wxy, S2.wxy = %e, %e\n", __func__, S1->wxy, S2->wxy);
  THEPRINTF("[%s]: S1.wyy, S2.wyy = %e, %e\n", __func__, S1->wyy, S2->wyy);
}

// upper borders excluded from "inside"
bool point_inside_box(double px, 
                      double py, 
                      const tPoint* cb, 
                      double hbw)
{
  const double cx = cb->x;
  const double cy = cb->y;
  return (px >= cx - hbw && px < cx + hbw && py >= cy - hbw && py < cy + hbw);
}

// theta = ( (box size) / (distance to box centre) ) ^ 2
double calc_theta(double px, 
                  double py, 
                  const tPoint* cb, 
                  double hbw)
{
  const double dx = cb->x - px;
  const double dy = cb->y - py;
  return (hbw * hbw / (dx * dx + dy * dy));
}

// Does box1 (center c1, half-width hw1) and box2 (c2, hw2) overlap in any way?
bool box_overlap_box(const tPoint* c1, 
                     double hw1, 
                     const tPoint* c2, 
                     double hw2) 
{
  const double right1 = c1->x + hw1;
  const double left1 = c1->x - hw1;
  const double right2 = c2->x + hw2;
  const double left2 = c2->x - hw2;
  if (right2 < left1) return false;
  if (right1 < left2) return false;
  const double top1 = c1->y + hw1;
  const double down1 = c1->y - hw1;
  const double top2 = c2->y + hw2;
  const double down2 = c2->y - hw2;
  if (top2 < down1) return false;
  if (top1 < down2) return false;
  return true;
}

// no error checking; internal call only
tValueTriad quadtree_leafsum_(const tQuadTree* root, 
                              const tPoint* target,
                              double epsq,
                              PotentialFuncPtr pfunc)
{
  tValueTriad sum = {0.0, 0.0, 0.0};
  for (int i = 0; i < root->npts; i++) {
    //const tValueTriad vt = logr_potential_(target->x, target->y, root->p[i].x, root->p[i].y, epsq);
    const tValueTriad vt = (*pfunc)(target->x, target->y, root->p[i].x, root->p[i].y, epsq);
    const double weight = root->p[i].w;
    sum.value += vt.value * weight;
    sum.gradx += vt.gradx * weight;
    sum.grady += vt.grady * weight;
  }
  return sum;
}

tValueTriad quadtree_sum_at_point(const tQuadTree* root,
                                  const tPoint* target,
                                  double theta,
                                  double epsq,
                                  PotentialFuncPtr pfunc,
                                  QuadrupoleFuncPtr qfunc)
{
  const bool target_inside_node = point_inside_box(target->x, 
                                                   target->y, 
                                                   &(root->cb), 
                                                   root->hbw);
  const bool isLeaf = (root->p != NULL); 

  if (target_inside_node && isLeaf) {
    // if the query point is inside the node and the node is a leaf: direct sum of the payload to target.
    return quadtree_leafsum_(root, target, epsq, pfunc);
  }

  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};

  if (target_inside_node && !isLeaf) { // always !isLeaf
    // if the query point is inside and the node is not a leaf: recursive calls with each child as new root.
    tValueTriad sum = {0.0, 0.0, 0.0};
    for (int c = 0; c < 4; c++) {
      if (child[c] == NULL) continue;
      const tValueTriad vtc = quadtree_sum_at_point(child[c], target, theta, epsq, pfunc, qfunc);
      sum.value += vtc.value;
      sum.gradx += vtc.gradx;
      sum.grady += vtc.grady;
    }
    return sum;
  }

  // Target point is outside this node
  const double theta_node = calc_theta(target->x, 
                                       target->y, 
                                       &(root->cb), 
                                       root->hbw);

  if (theta_node < theta) {
    // if the query point is outside the node: if "far-field" (w.r.t. theta) then use quadrupole approximation
    const double X = target->x - root->cb.x;
    const double Y = target->y - root->cb.y;
    const tValueTriad vt = (*qfunc)(X, Y, &(root->stats));
    return vt;
  }

  // target point is not "far-field" but it is outside this node; could be a leaf
  if (isLeaf) {
    return quadtree_leafsum_(root, target, epsq, pfunc);
  }

  // not far field and not leaf; recurse
  tValueTriad sum = {0.0, 0.0, 0.0};
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL) continue;
    const tValueTriad vtc = quadtree_sum_at_point(child[c], target, theta, epsq, pfunc, qfunc);
    sum.value += vtc.value;
    sum.gradx += vtc.gradx;
    sum.grady += vtc.grady;
  }

  return sum;
}

// 00 = 3rd quadrant (mm)
// 01 = 2nd quadrant (pm)
// 10 = 4th quadrant (mp)
// 11 = 1st quadrant (pp)
uint8_t quadrant_index(const tPointPayload* p, 
                       const tPoint* c) 
{
  uint8_t k = 0;
  if (p->x >= c->x) k |= 0x01;
  if (p->y >= c->y) k |= 0x02;
  return k;
}

// return the number of nodes needed during the construction process (does not include root)
// point is the set of points (it is rearranged internally & accessed through tree traversal later)
// scratch is a temporary sorting area
int build_quadtree(tQuadTree* root, 
                   int maxInLeaf,
                   int maxDepth,
                   int level,
                   int np,
                   tPointPayload* point, 
                   tPointPayload* scratch,
                   int nn, 
                   tQuadTree* node) 
{
  // On entry, assume that root->cb and root->hbw have been set already... and go from there

  root->pp = NULL;
  root->pm = NULL;
  root->mm = NULL;
  root->mp = NULL;

  // Are we a leaf node? If so, then mark node; compute its stats w.r.t. to node center; and stop
  if (np <= maxInLeaf || level == maxDepth) {
    root->npts = np;
    root->p = point;
    calc_source_summary(&(root->stats), root->npts, root->p, root->cb.x, root->cb.y);
    //THEPRINTF("[%s]: cbx, cby, hbw = %f, %f, %f (leaf, np = %i)\n", __func__, root->cb.x, root->cb.y, root->hbw, root->npts);
    return 0;
  }

  // Too many points to store as a leaf; subdivide by recursion..
  // First apply counting sort on given point array and then 
  // call quadtree builder on each required (non-empty) quadrant point set.

  root->npts = 0;
  root->p = NULL;

  if (nn == 0)
    return 0; // cannot continue since there are no more nodes to grab from the heap

  int count[4] = {0, 0, 0, 0}; // mm, pm, mp, pp

  for (int i = 0; i < np; i++) {
    const uint8_t ikey = quadrant_index(&(point[i]), &(root->cb));
    count[ikey]++;
    scratch[i] = point[i];
  }

  for (int i = 1; i <= 3; i++) {
    count[i] += count[i - 1]; // will hold offsets in the end
  }

  for (int i = np - 1; i >= 0; i--) {
    const uint8_t ikey = quadrant_index(&(scratch[i]), &(root->cb));
    count[ikey]--;
    point[count[ikey]] = scratch[i];
  }

  // point array re-ordered!
  // index: count[0] .. count[1] - 1 : 3rd
  // index: count[1] .. count[2] - 1 : 2nd
  // index: count[2] .. count[3] - 1 : 4th
  // index: count[3] .. np - 1       : 1st

  int n_3rd = count[1] - count[0];
  int n_2nd = count[2] - count[1];
  int n_4th = count[3] - count[2];
  int n_1st = np - count[3];

  // OK now for each quadrant set which is non-empty; create a child node and
  // recurse; making sure the available node store is updated ..

  // reset the stats for the non-leaf nodes; aggregate child stats below
  memset(&(root->stats), 0 , sizeof(tMomentStats));

  /*THEPRINTF("[%s]: pp, pm, mm, mp = %i, %i, %i, %i (sum = %i)\n", 
    __func__, n_1st, n_2nd, n_3rd, n_4th, n_1st + n_2nd + n_3rd + n_4th); */

  const double hhbw = root->hbw / 2.0;

  int nnodes = 0;

  if (n_3rd != 0) {
    nnodes++;
    root->mm = node++;
    nn--;
    root->mm->cb.x = root->cb.x - hhbw;
    root->mm->cb.y = root->cb.y - hhbw;
    root->mm->hbw = hhbw;
    int nn_3rd = build_quadtree(root->mm, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_3rd,
                                &(point[count[0]]), 
                                &(scratch[count[0]]),
                                nn, 
                                node);
    nnodes += nn_3rd;
    nn -= nn_3rd;
    node += nn_3rd;

    add_source_translation(&(root->stats), 
                           &(root->mm->stats),
                           -hhbw,
                           -hhbw);
  }

  if (n_2nd != 0) {
    nnodes++;
    root->pm = node++;
    nn--;
    root->pm->cb.x = root->cb.x + hhbw;
    root->pm->cb.y = root->cb.y - hhbw;
    root->pm->hbw = hhbw;
    int nn_2nd = build_quadtree(root->pm, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_2nd,
                                &(point[count[1]]), 
                                &(scratch[count[1]]),
                                nn, 
                                node);
    nnodes += nn_2nd;
    nn -= nn_2nd;
    node += nn_2nd;

    add_source_translation(&(root->stats), 
                           &(root->pm->stats),
                           hhbw,
                           -hhbw);
  }

  if (n_4th != 0) {
    nnodes++;
    root->mp = node++;
    nn--;
    root->mp->cb.x = root->cb.x - hhbw;
    root->mp->cb.y = root->cb.y + hhbw;
    root->mp->hbw = hhbw;
    int nn_4th = build_quadtree(root->mp, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_4th,
                                &(point[count[2]]), 
                                &(scratch[count[2]]),
                                nn, 
                                node);
    nnodes += nn_4th;
    nn -= nn_4th;
    node += nn_4th;

    add_source_translation(&(root->stats), 
                           &(root->mp->stats),
                           -hhbw,
                           hhbw);
  }

  if (n_1st != 0) {
    nnodes++;
    root->pp = node++;
    nn--;
    root->pp->cb.x = root->cb.x + hhbw;
    root->pp->cb.y = root->cb.y + hhbw;
    root->pp->hbw = hhbw;
    int nn_1st = build_quadtree(root->pp, 
                                maxInLeaf,
                                maxDepth,
                                level + 1,
                                n_1st,
                                &(point[count[3]]), 
                                &(scratch[count[3]]),
                                nn, 
                                node);
    nnodes += nn_1st;
    nn -= nn_1st;
    node += nn_1st;

    add_source_translation(&(root->stats), 
                           &(root->pp->stats),
                           hhbw,
                           hhbw);
  }

  return nnodes;
}

// Simple box-query that just returns the total number of points found
int quadtree_box_query_count(const tQuadTree* root, 
                             const tPoint* cq, 
                             double hwq)
{
  if (!box_overlap_box(&(root->cb), root->hbw, cq, hwq))
    return 0;

  if (root->p != NULL) {
    int nqpts = 0;
    for (int j = 0; j < root->npts; j++) {
      if (point_inside_box(root->p[j].x, root->p[j].y, cq, hwq)) {
        nqpts++;
      }
    }
    return nqpts;
  }

  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int ncq = 0;

  for (int c = 0; c < 4; c++) {
    if (child[c] != NULL)
      ncq += quadtree_box_query_count(child[c], cq, hwq);
  }

  return ncq;
}

// Count the number of nodes in the tree (supposed to match the return value of the build code)
int count_quadtree_nodes(const tQuadTree* root) {
  if (root->p != NULL)
    return 1;  // this is a leaf node
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int nn = 1;  // count this node.. plus all the nodes within the children..
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    nn += count_quadtree_nodes(child[c]);
  }
  return nn;
}

// Tally up the total number of points stored in the leaves.
int count_quadtree_points(const tQuadTree* root) {
  if (root->p != NULL)
    return root->npts;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int n = 0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    n += count_quadtree_points(child[c]);
  }
  return n;
}

// Return total number of leaf nodes
int count_quadtree_leaves(const tQuadTree* root) {
  if (root->p != NULL)
    return 1;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int n = 0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    n += count_quadtree_leaves(child[c]);
  }
  return n;
}

// Tally up maximum depth of tree.
int count_maximum_depth(const tQuadTree* root, int ref) {
  if (root->p != NULL)
    return ref;
  int d = ref;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    const int thisd = count_maximum_depth(child[c], ref + 1);
    if (thisd > d)
      d = thisd;
  }
  return d;
}

// tally up the average depth of the tree
double count_average_depth(const tQuadTree* root, double ref) {
  if (root->p != NULL)
    return ref;
  const tQuadTree* child[4] = {root->pp, root->pm, root->mm, root->mp};
  int nc = 0;
  double dsum = 0.0;
  for (int c = 0; c < 4; c++) {
    if (child[c] == NULL)
      continue;
    const double thisd = count_average_depth(child[c], ref + 1.0);
    dsum += thisd;
    nc++;
  }
  return (dsum / nc);
}

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
#define ARG_EP       prhs[3]
#define ARG_MAXLEAF  prhs[4]
#define ARG_EPK      prhs[5]
#define ARG_POTSTR   prhs[6]

#define MAX_OUTARGS  2
#define OUT_WP       plhs[0]
#define OUT_NNODES   plhs[1]

void mexFunction(int nlhs, 
                 mxArray** plhs,
                 int nrhs, 
                 const mxArray** prhs)
{
  if (nrhs != NUM_INARGS) {
    THEPRINTF("USAGE: wp = qtreehat(xp, yp, vp, ep, maxleaf, epk, pot);\n");
    THEERRMSG("Incorrect number of input arguments provided");
  }

  if (nlhs > MAX_OUTARGS) {
    THEERRMSG("Too many output arguments requested.");
  }

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

  if (!(isDoubleRealVector(ARG_EP) && mxGetNumberOfElements(ARG_EP) == 1)) {
    THEERRMSG("ep must be a real-valued scalar"); 
  }

  const double EP = *mxGetPr(ARG_EP);
  const bool EP_is_valid = (EP > 0.0 & EP < 1.0);

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
    THEERRMSG("pot must be either \"log\" or \"inv\" (case insensitive)");
  }

  const double epksq = EPK * EPK;

  PotentialFuncPtr pfunc  = (use_log_pot ? &logr_potential_ : &invr_potential_);
  QuadrupoleFuncPtr qfunc = (use_log_pot ? &eval_logr_quadrupole_ : &eval_invr_quadrupole_);

  const int Np = nxp;
  OUT_WP = mxCreateDoubleMatrix(Np, 3, mxREAL);
  double* pWP = mxGetPr(OUT_WP);

  if (!EP_is_valid) {
    for (int i = 0; i < Np; i++) {
      const double xi = pXP[i];
      const double yi = pYP[i];
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
      pWP[i + Np] = sumfx;
      pWP[i + 2 * Np] = sumfy;
    }

    return;
  }

  const int maxNumberOfNodes = Np * 4; // (!!! !!! !!!) probably only 2 needed here

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
    THEPRINTF("[%s]: nno = %i\n", __func__, nno);

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
  for (int i = 0; i < Np; i++) {
    const double xi = pXP[i];
    const double yi = pYP[i];
    const tPoint ith_query = {xi, yi};
    const tValueTriad wi = quadtree_sum_at_point(&rootNode, 
                                                 &ith_query, 
                                                 theta, 
                                                 epksq, 
                                                 pfunc, 
                                                 qfunc);
    pWP[i] = wi.value;
    pWP[i + Np] = wi.gradx;
    pWP[i + 2 * Np] = wi.grady;
  }

  THEFREE(pt);
  THEFREE(pt_scratch);
  THEFREE(nodes);

  if (nlhs == 2) {
    OUT_NNODES = mxCreateDoubleScalar((double)(nno + 1));
  }

  return;
}

// TODO: optional set of target points (otherwise evaluate at source points)
// TODO: switchable level of detail: 0th, 1st, 2nd? i.e. stop at dipole or go full quadrupole?
// TODO: internally clock the non-allocation steps in the calculations.. return in a 2nd output
// TODO: review implementation of logr quadrupole (use Xhat, Yhat version?)
// TODO: what is the maximum possible number of nodes in the quadtree as function of maxleaf and numpoints?
