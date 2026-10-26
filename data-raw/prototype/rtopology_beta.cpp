// Prototype (not part of the package): rtopology_beta() in C++ through Rcpp,
// written to reproduce the R version exactly for the same seed. It uses R's
// own RNG (unif_rand), lgammafn and revsort, and the steps of R's
// sample(x, 1, prob = p) without replacement (FixupProb, then
// ProbSampleNoReplace), so that the same uniforms map to the same splits.
#include <Rcpp.h>
#include <R_ext/Utils.h>   // revsort
#include <Rmath.h>         // lgammafn
#include <vector>
using namespace Rcpp;

// one draw of sample(seq_len(len), 1, prob = w), as R does it
static int sample_one(std::vector<double>& w, std::vector<int>& perm, int len) {
  double sum = 0.0;
  for (int i = 0; i < len; i++) sum += w[i];           // FixupProb
  for (int i = 0; i < len; i++) w[i] /= sum;
  for (int i = 0; i < len; i++) perm[i] = i + 1;       // ProbSampleNoReplace
  revsort(w.data(), perm.data(), len);
  double rT = unif_rand(), mass = 0.0;
  int j;
  for (j = 0; j < len - 1; j++) {
    mass += w[j];
    if (rT <= mass) break;
  }
  return perm[j];
}

// [[Rcpp::export]]
IntegerMatrix rtopology_edges_cpp(int n, double beta) {
  IntegerMatrix edge(2 * n - 2, 2);
  int n_edge = 0, next_node = n + 1;
  std::vector<int> st_start{1}, st_m{n}, st_parent{0};
  std::vector<double> lw(n), w(n);
  std::vector<int> perm(n);
  GetRNGstate();
  while (!st_start.empty()) {
    int start = st_start.back(), m = st_m.back(), parent = st_parent.back();
    st_start.pop_back(); st_m.pop_back(); st_parent.pop_back();
    int id;
    if (m == 1) id = start; else id = next_node++;
    if (parent > 0) { edge(n_edge, 0) = parent; edge(n_edge, 1) = id; n_edge++; }
    if (m > 1) {
      int len = m - 1;
      double mx = R_NegInf;
      for (int i = 1; i <= len; i++) {
        lw[i - 1] = R::lgammafn(beta + i + 1) + R::lgammafn(beta + m - i + 1) -
                    R::lgammafn(i + 1) - R::lgammafn(m - i + 1);
        if (lw[i - 1] > mx) mx = lw[i - 1];
      }
      for (int i = 0; i < len; i++) w[i] = exp(lw[i] - mx);
      int k = sample_one(w, perm, len);
      st_start.push_back(start + k); st_m.push_back(m - k); st_parent.push_back(id);
      st_start.push_back(start);     st_m.push_back(k);     st_parent.push_back(id);
    }
  }
  PutRNGstate();
  return edge;
}
