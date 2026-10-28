#include <algorithm>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "../include/mvpolynomial.hpp"

/*
 * SCIPpp and SCIP-SDP dependencies
 */
#include <scippp/model.hpp>                 /* scippp::Model, Var, LinExpr, Solution */
#include <scip/cons_sdp.h>                  /* SCIP-SDP: SCIPcreateConsSdp          */
#include <scip/scipsdpdefplugins.h>         /* SCIP-SDP: SCIPSDPincludeDefaultPlugins */

/*
 * Last convergent h/k of the continued fraction of x with k <= maxden.
 * Exact integer arithmetic (Euclid on num/den): no floating point error
 * and no overflow, for any size of x.
 *   h_n = a_n h_{n-1} + h_{n-2},   k_n = a_n k_{n-1} + k_{n-2}
 */
static mpq_class 
cf_approx(const mpq_class& x, const mpz_class& maxden) {
  mpz_class n = x.get_num(), d = x.get_den();
  bool neg = (n < 0);
  if (neg) n = -n;

  mpz_class hm2 = 0, hm1 = 1;   /* h_{-2}, h_{-1} */
  mpz_class km2 = 1, km1 = 0;   /* k_{-2}, k_{-1} */
  mpz_class a, h, k, r;

  while (d != 0) {
    a = n / d;                  /* floor, since n >= 0 and d > 0 */
    h = a * hm1 + hm2;
    k = a * km1 + km2;
    if (k > maxden) break;
    hm2 = hm1; hm1 = h;
    km2 = km1; km1 = k;
    r = n - a * d;
    n = d;
    d = r;
  }
  mpq_class res(hm1, km1);
  res.canonicalize();
  return neg ? mpq_class(-res) : res;
}

/* restores the default mpf precision when leaving refine() */
struct MpfPrecGuard {
  mp_bitcnt_t old;
  explicit MpfPrecGuard(mp_bitcnt_t bits) : old(mpf_get_default_prec()) { mpf_set_default_prec(bits); }
  ~MpfPrecGuard() { mpf_set_default_prec(old); }
};

/* "0.123e-45" style formatting for mpf (safe below the double range) */
static std::string 
sci(const mpf_class& x) {
  if (x == 0) return std::string("0");
  mp_exp_t e;
  std::string s = x.get_str(e, 10, 3);
  return "0." + s + "e" + std::to_string(static_cast<long>(e));
}

/* Gauss elimination with partial pivoting; solves G x = b (result in b) */
static bool 
solve_linear(std::vector<std::vector<mpf_class>>& G, std::vector<mpf_class>& b,
                  const mpf_class& ptol) {
  const size_t m = b.size();
  for (size_t c = 0; c < m; c++) {
    size_t piv = c;
    mpf_class best = abs(G[c][c]);
    for (size_t r = c + 1; r < m; r++) {
      mpf_class v = abs(G[r][c]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < ptol) return false;
    if (piv != c) { std::swap(G[piv], G[c]); std::swap(b[piv], b[c]); }
    for (size_t r = c + 1; r < m; r++) {
      if (G[r][c] == 0) continue;
      mpf_class f = G[r][c] / G[c][c];
      for (size_t k = c; k < m; k++) G[r][k] -= f * G[c][k];
      b[r] -= f * b[c];
    }
  }
  for (size_t c = m; c-- > 0; ) {
    mpf_class s = b[c];
    for (size_t k = c + 1; k < m; k++) s -= G[c][k] * b[k];
    b[c] = s / G[c][c];
  }
  return true;
}

mvpolyT::mvpolyT () {
  std::cin>>nvars;
  double temp_c;
  exps e(nvars);
  ismbcomputed = false;
  verbose = true;
  
  while(std::cin>>temp_c) {
    coeff c(temp_c);
    
    c.canonicalize();
    for (unsigned int i = 0; i<nvars; i++)
      std::cin>>e[i];
    addterm(e,c);
  }
}

void 
mvpolyT::addterm (const exps& _e, coeff _c) {
  if(_c == 0)
    return;

  auto it = terms.find(_e);
  if(it != terms.end()) {
    it->second += _c;

    if(it->second == 0)
      terms.erase(it);
  }
  else {
    terms[_e] = _c; 
  }
}

unsigned int 
mvpolyT::degree() const {
  int d = 0;
  
  std::map<exps, coeff>::const_iterator it;
  for (it = terms.begin(); it!= terms.end(); it++) {
    int sum = 0;
    for (unsigned int e : it->first)
      sum+=e;

    d = std::max(d,sum);
  }
  return(d);
}

bool 
mvpolyT::homogeneous() const {
  if(terms.empty())
    return (true);

  std::map<exps, coeff>::const_iterator it;
  it = terms.begin();
  unsigned int firs_degree = 0;
  for (unsigned int e: it->first)
    firs_degree+=e;

  for (it = terms.begin(); it != terms.end(); it++) {
    unsigned int d = 0;
    for (unsigned int e: it->first)
      d+=e;

    if(d != firs_degree)
      return(false);
  }
  return(true);
}

/*
 * The number of entries in this vector is C(n+d, d).
 * Note that the monomial indices are all tuples of 
 * length n with non-negative integrer entries summing
 * to at most d (which can be generated with recursion
 * and backtracking).
 * Remark: If the input polynomial is homogeneous of 
 * degree 2d, then it is sufficiente to restrict the 
 * components of monomialbasis to the monomials of 
 * degree exactly equal to d.
 */ 
void 
mvpolyT::rec_monomialbasis(unsigned int v, unsigned int remain, exps& curr, std::vector<exps>& ans) const {
  /*
   * Last variable
   */
  if (v == nvars-1) {
    curr[v] = remain;
    ans.push_back(curr);
    return;
  }

  /*
   * For the actual variable, try all the powers
   */
  for (unsigned int e = 0; e<=remain; e++) {
    curr[v]=e;
    rec_monomialbasis(v+1, remain - e, curr, ans);
  }
}

// TODO Newton politope optimization
const monomialbasisT& 
mvpolyT::monomialbasis() {
  if(!ismbcomputed) {
    unsigned int twod = degree();
    if (twod % 2 != 0)
      throw std::runtime_error("Polynomials with odd degree always have negative points");
  
    mb.degree = twod / 2;
    mb.homogeneous = homogeneous();

    if (mb.homogeneous) {
      exps curr(nvars);
      rec_monomialbasis(0, mb.degree, curr, mb.monomials);
    }
    else {
      for (unsigned int d = 0; d<= mb.degree; d++) {
	    exps curr(nvars);
	    rec_monomialbasis(0, d, curr, mb.monomials);
      }
    }
  
    mb.size = mb.monomials.size();
    ismbcomputed = true;
    }
  return(mb);
}

/*
 * Given a monomial basis, does not matter the ordering this function
 * uses the gram's matrix method for computing a vector of linear
 */
std::vector<constraintT> 
mvpolyT::gram() {
  const monomialbasisT& mb = monomialbasis();
  std::vector<constraintT> ans;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> entryindex; 
  /* 
   * Each (i,j) mapsto a decision variable say y_k, with index k.
   * We have N(N+1)/2 decision variables since the Gram matrix is symetric.
  */
  unsigned int index = 0;
  for(unsigned int i = 0; i<mb.size; i++) {
    for(unsigned int j = i; j<mb.size; j++) {
      index++;
      entryindex[{i,j}] = index;
    }
  }

  std::map<exps, std::vector<std::pair<unsigned int, coeff>>> monomialindex; 
  /* 
   * Each variable needs to be compared to a coefficient and it is 
   * dicted with the monomial. Expoents maps to the index of decision
   * variable plus the coefficient, note that since the Gram matrix is symetric,
   * coeff = 1 if i==j and coeff = 2 if i != j.
   */
  for(unsigned int i = 0; i<mb.size; i++) {
    for(unsigned int j = i; j<mb.size; j++) {
      exps product(nvars);
      for (unsigned int k = 0; k < nvars; k++) {
        product[k] = mb.monomials[i][k] + mb.monomials[j][k];
      }

      unsigned int y_k = entryindex[{i,j}];
      coeff c = (i==j)? 1 : 2;

      monomialindex[product].push_back({y_k, c});
    }
  }

  // Generating linear constraints
  for (auto const& pair : monomialindex) {
    constraintT ctr;
    auto term_it = terms.find(pair.first);
    ctr.rhs = (term_it != terms.end()) ? term_it->second : coeff(0);

    for(const auto& ic : pair.second)
      ctr.lhs[ic.first] = ic.second;

    ans.push_back(ctr);
  }

  return(ans);
}

/*
 * Builds and solves the Gram-matrix feasibility SDP in memory:
 */
void
mvpolyT::solvesdp() {
  const monomialbasisT& mb = monomialbasis();
  const std::vector<constraintT> ctrs = gram();
  const size_t N  = mb.size;
  const size_t ny = N * (N + 1) / 2;            /* distinct Gram entries */

  SCIP* scip = nullptr;
  if (SCIPcreate(&scip) != SCIP_OKAY)
    throw std::runtime_error("SCIPcreate failed");
  if (SCIPSDPincludeDefaultPlugins(scip) != SCIP_OKAY)
    throw std::runtime_error("SCIPSDPincludeDefaultPlugins failed");

  bool optimal = false;

  { /* scope: the Model must be destroyed before SCIPfree(&scip) */
    scippp::Model model("gram", scip, /*includeDefaultPlugins=*/false);

    /*Vector of variables, it is diferent of the vector from SCIPpp which needs a know size*/
    std::vector<scippp::Var> y;
    y.reserve(ny);
    for (size_t k = 0; k < ny; k++)
      y.push_back(model.addVar("y" + std::to_string(k), 0.0, scippp::VarType::CONTINUOUS,
                               std::nullopt, std::nullopt));

    for (size_t c = 0; c < ctrs.size(); c++) {
      scippp::LinExpr e;
      for (const auto& [idx1, a] : ctrs[c].lhs)
        e += a.get_d() * y[idx1 - 1];
      model.addConstr(e == ctrs[c].rhs.get_d(), "gram_" + std::to_string(c));
    }

    /*
     * SCIP-SDP integration, sdp constraints generation
     */
    std::vector<int> rows(ny), cols(ny), nvarnonz(ny, 1);   /* 1 nonzero per variable */
    std::vector<double> vals(ny, 1.0);
    std::vector<int*> rowp(ny), colp(ny);
    std::vector<double*> valp(ny);
    std::vector<SCIP_VAR*> scipvars(ny);

    size_t k = 0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i; j < N; j++) {
        rows[k] = static_cast<int>(j);
        cols[k] = static_cast<int>(i);
        rowp[k] = &rows[k];
        colp[k] = &cols[k];
        valp[k] = &vals[k];
        scipvars[k] = y[k].getVar();
        k++;
      }
    }

    SCIP_CONS* sdpcons = nullptr;
    if (SCIPcreateConsSdp(scip, &sdpcons, "gram_psd",
                          static_cast<int>(ny), static_cast<int>(ny), static_cast<int>(N),
                          nvarnonz.data(), colp.data(), rowp.data(), valp.data(),
                          scipvars.data(),
                          0, nullptr, nullptr, nullptr,
                          FALSE) != SCIP_OKAY)
      throw std::runtime_error("SCIPcreateConsSdp failed");
    if (SCIPaddCons(scip, sdpcons) != SCIP_OKAY)
      throw std::runtime_error("SCIPaddCons failed");
    SCIPreleaseCons(scip, &sdpcons);

    model.solve();
    optimal = (model.getStatus() == SCIP_STATUS_OPTIMAL && model.getNSols() > 0);

    if (optimal) {
      scippp::Solution sol = model.getBestSol();
      solution.y.assign(ny, 0.0);
      for (size_t kk = 0; kk < ny; kk++)
        solution.y[kk] = y[kk].getSolVal(sol);
    }
  } /* Model destroyed here */

  SCIPfree(&scip);

  if (!optimal)
    throw std::runtime_error("SCIP-SDP found no optimal solution");

  solution.y_exact.clear();
  solution.y_refined.clear();
  solution.refined = false;
  solution.nvars  = static_cast<unsigned int>(ny);
  solution.size   = N;
  solution.solved = true;
}

/*
 * Rationalization with continued fractions: every y_k (the refined
 * mpf value if available, else the double from the solver) is replaced
 * by its last convergent with denominator <= maxden. 
 */
void
mvpolyT::continued_fractions(const mpz_class& maxden) {
  const size_t n = solution.y.size();
  const bool useref = solution.refined && solution.y_refined.size() == n;
  solution.y_exact.assign(n, coeff(0));

  for (size_t k = 0; k < n; k++) {
    mpq_class x;
    if (useref)
      mpq_set_f(x.get_mpq_t(), solution.y_refined[k].get_mpf_t()); /* exact */
    else
      x = mpq_class(solution.y[k]);                                /* exact */
    x.canonicalize();
    solution.y_exact[k] = cf_approx(x, maxden);
  }
}

/*
 * Projection into the affine subspace of valid gram matrices
 * (Peyrl-Parrilo, Prop. 7). Entries of the same monomial gamma form disjoint
 * groups, and each y_k in a group is shifted by e_gamma / n(gamma)
 * where n(gamma) = sum of the multiplicities (1 diagonal / 2 off-diagonal).
 */
// TODO data structure for monomial hash table
void
mvpolyT::projection() {
  if(solution.y_exact.empty())
    throw std::runtime_error("It is necessary to solve first");

  const monomialbasisT& mb = monomialbasis();

  /*
   * expoents |-> y index and coefficient
   */
  std::map<exps, std::vector<std::pair<unsigned int, coeff>>> monomialindex;
  unsigned int k = 0;
  for (unsigned int i = 0; i < mb.size; i++) {
    for (unsigned int j = i; j < mb.size; j++) {
      exps product(nvars);
      for (unsigned int v = 0; v < nvars; v++)
        product[v] = mb.monomials[i][v] + mb.monomials[j][v];
      coeff c = (i == j) ? 1 : 2;
      monomialindex[product].push_back({k, c});
      k++;
    }
  }

  coeff maxerr = 0;
  for (const auto& item : monomialindex) {
    const exps& gamma = item.first;
    const auto& entries = item.second;

    coeff n_gamma = 0;
    coeff current_val = 0;
    for (const auto& entry : entries) {
      n_gamma += entry.second;
      current_val += entry.second * solution.y_exact[entry.first];
    }

    coeff target_val = 0;
    auto target_it = terms.find(gamma);
    if (target_it != terms.end())
      target_val = target_it->second;

    coeff error = current_val - target_val;
    if (abs(error) > maxerr) maxerr = abs(error);

    if (error != 0) {
      coeff correction = error / n_gamma;
      for (const auto& entry : entries) {
        solution.y_exact[entry.first] -= correction;
        solution.y_exact[entry.first].canonicalize();
      }
    }
  }
  if (verbose)
    std::cout << "[projection] max coefficient error before projection = "
              << maxerr.get_d() << '\n';
}

/*
 * Gauss-Newton refinement.
 *
 * Write the numerical Gram matrix as Q = A^T A (Cholesky), i.e. f ~ sum_i p_i^2
 * with p_i = sum_j A[i][j] z_j. Look for dA with
 *        f - sum_i (p_i + dp_i)^2 = 0   linearized:   2 sum_i p_i dp_i = r,
 *        r = f - sum_i p_i^2.
 * It is an underdetermined linear system J dA = r (J has m rows = #monomials
 * of degree <= 2d); the minimum norm solution is dA = J^T (J J^T)^{-1} r.
 * The new residual is -sum dp_i^2 = O(|dA|^2): quadratic convergence.
 *
 * W = A^T A is positive semidefinite by construction and f - z^T W z -> 0,
 * so the refined matrix is (almost) exactly in the affine subspace while
 * staying inside the PSD cone. That is what makes the later projection safe
 * (delta -> 0 in Peyrl-Parrilo's Prop. 8).
 *
 * Assumes Q is positive definite (strict feasibility), i.e. rank N.
 * The refined matrix is stored in solution.y_refined only if it is better
 * than the solver's one (solution.refined tells whether it was stored).
 */
void
mvpolyT::refine(unsigned int prec_bits, unsigned int tol_digits, unsigned int maxiter) {
  if(!solution.solved)
    throw std::runtime_error("It is necessary to solve first");

  const monomialbasisT& B = monomialbasis();
  const size_t N = B.size;
  const size_t ny = N * (N + 1) / 2;
  MpfPrecGuard guard(prec_bits);

  /* index of every monomial gamma = m_j + m_l */
  std::map<exps, size_t> gidx;
  std::vector<std::vector<size_t>> pg(N, std::vector<size_t>(N));
  for (size_t j = 0; j < N; j++) {
    for (size_t l = 0; l < N; l++) {
      exps g(nvars);
      for (unsigned int v = 0; v < nvars; v++)
        g[v] = B.monomials[j][v] + B.monomials[l][v];
      size_t id = gidx.size();
      pg[j][l] = gidx.emplace(g, id).first->second;
    }
  }
  const size_t m = gidx.size();

  std::vector<mpf_class> target(m);
  for (const auto& t : terms) {
    auto it = gidx.find(t.first);
    if (it == gidx.end())
      continue;
    mpf_set_q(target[it->second].get_mpf_t(), t.second.get_mpq_t());
  }

  /* numerical Gram matrix */
  std::vector<std::vector<mpf_class>> Q(N, std::vector<mpf_class>(N));
  {
    size_t k = 0;
    for (size_t i = 0; i < N; i++)
      for (size_t j = i; j < N; j++, k++)
        Q[i][j] = Q[j][i] = solution.y[k];
  }

  /* Q + eta*I = A^T A, A upper triangular (Cholesky). eta > 0 only if the solver's Q
     is numerically indefinite; Gauss-Newton absorbs that perturbation. */
  std::vector<std::vector<mpf_class>> A(N, std::vector<mpf_class>(N));
  mpf_class eta = 0;
  for (int attempt = 0; ; attempt++) {
    bool ok = true;
    for (size_t i = 0; i < N; i++) for (size_t j = 0; j < N; j++) A[i][j] = 0;
    for (size_t i = 0; i < N && ok; i++) {
      mpf_class s = Q[i][i] + eta;
      for (size_t kk = 0; kk < i; kk++) s -= A[kk][i] * A[kk][i];
      if (s <= 1e-13 * abs(Q[i][i])) { ok = false; break; }
      A[i][i] = sqrt(s);
      for (size_t j = i + 1; j < N; j++) {
        mpf_class t = Q[i][j];
        for (size_t kk = 0; kk < i; kk++) t -= A[kk][i] * A[kk][j];
        A[i][j] = t / A[i][i];
      }
    }
    if (ok) break;
    if (attempt >= 40)
      throw std::runtime_error("refine: Gram matrix is not positive definite (not strictly feasible?)");
    eta = (eta == 0) ? mpf_class(1e-12) : eta * 10;
  }
  if (verbose && eta != 0)
    std::cout << "[refine] diagonal shift eta = " << sci(eta) << " was needed\n";

  const mpf_class tol(std::string("1e-") + std::to_string(tol_digits));
  const mpf_class ptol(std::string("1e-") + std::to_string(prec_bits * 3 / 10 - 10));

  std::vector<std::vector<mpf_class>> Qa(N, std::vector<mpf_class>(N));
  std::vector<mpf_class> res(m);
  mpf_class theta0 = 0, theta = 0;
  bool converged = false;

  for (unsigned int it = 0; ; it++) {

    /* Qa = A^T A and residual r = f - sum_gamma (...) */
    for (size_t j = 0; j < N; j++) {
      for (size_t l = j; l < N; l++) {
        mpf_class s = 0;
        for (size_t i = 0; i < N; i++) s += A[i][j] * A[i][l];
        Qa[j][l] = Qa[l][j] = s;
      }
    }
    std::vector<mpf_class> c(m);
    for (size_t j = 0; j < N; j++)
      for (size_t l = 0; l < N; l++)
        c[pg[j][l]] += Qa[j][l];

    theta = 0;
    for (size_t g = 0; g < m; g++) {
      res[g] = target[g] - c[g];
      mpf_class v = abs(res[g]);
      if (v > theta) theta = v;
    }
    if (it == 0) theta0 = theta;
    if (verbose)
      std::cout << "[refine] iter " << it << "  backward error theta = " << sci(theta) << '\n';

    if (theta < tol) { converged = true; break; }
    if (it >= maxiter) break;

    /*
     * J J^T:  J[gamma][(i,j)] = 2 A[i][l] with m_j + m_l = gamma
     * (for fixed j, distinct l give distinct gamma)
     */
    std::vector<std::vector<mpf_class>> G(m, std::vector<mpf_class>(m));
    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++)
        for (size_t l1 = 0; l1 < N; l1++) {
          mpf_class a1 = 2 * A[i][l1];
          for (size_t l2 = 0; l2 < N; l2++)
            G[pg[j][l1]][pg[j][l2]] += a1 * (2 * A[i][l2]);
        }

    std::vector<mpf_class> lambda = res;
    if (!solve_linear(G, lambda, ptol)) {
      if (verbose) std::cout << "[refine] J J^T singular, stopping\n";
      break;
    }

    /* dA = J^T lambda */
    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++) {
        mpf_class s = 0;
        for (size_t l = 0; l < N; l++) s += A[i][l] * lambda[pg[j][l]];
        A[i][j] += 2 * s;
      }
  }

  if (!converged && verbose)
    std::cout << "[refine] warning: tolerance not reached\n";

  /* only keep the refined matrix if it is really better than the solver's one */
  if (theta < theta0 || converged) {
    solution.y_refined.assign(ny, mpf_class(0));
    size_t k = 0;
    for (size_t i = 0; i < N; i++)
      for (size_t j = i; j < N; j++, k++)
        solution.y_refined[k] = Qa[i][j];
    solution.refined = true;
  }
}

/*
 * Exact identity  f == z^T Q z  with Q = y_exact (checks every coefficient).
 */
bool
mvpolyT::verify_identity() const {
  std::map<exps, coeff> acc;
  size_t k = 0;
  for (unsigned int i = 0; i < mb.size; i++) {
    for (unsigned int j = i; j < mb.size; j++, k++) {
      exps p(nvars);
      for (unsigned int v = 0; v < nvars; v++)
        p[v] = mb.monomials[i][v] + mb.monomials[j][v];
      acc[p] += ((i == j) ? 1 : 2) * solution.y_exact[k];
    }
  }
  for (const auto& a : acc) {
    auto it = terms.find(a.first);
    coeff t = (it == terms.end()) ? coeff(0) : it->second;
    if (a.second != t) return false;
  }
  for (const auto& t : terms)
    if (acc.find(t.first) == acc.end()) return false;
  return true;
}

std::vector<std::vector<coeff>>
mvpolyT::gram_exact() const {
  std::vector<std::vector<coeff>> Q(mb.size, std::vector<coeff>(mb.size, 0));
  size_t k = 0;
  for (unsigned int i = 0; i < mb.size; i++)
    for (unsigned int j = i; j < mb.size; j++, k++)
      Q[i][j] = Q[j][i] = solution.y_exact[k];
  return Q;
}

/*
 * Exact LDL^T without pivoting. Returns whether Q is PSD:
 *   - every D[i] >= 0, and
 *   - if D[i] == 0 the whole remaining column must be 0
 *     (otherwise Q is not PSD).
 */
bool
mvpolyT::LDLT_factor(const std::vector<std::vector<coeff>>& Q,
                     std::vector<std::vector<coeff>>& L,
                     std::vector<coeff>& D) const {
  const size_t n = Q.size();
  L.assign(n, std::vector<coeff>(n, 0));
  D.assign(n, 0);
  bool psd = true;

  for (size_t i = 0; i < n; i++)
    L[i][i] = 1;
  for (size_t i = 0; i < n; i++) {
    coeff sum = 0;
    for (size_t j = 0; j < i; j++)
      sum += L[i][j] * L[i][j] * D[j];
    D[i] = Q[i][i] - sum;
    if (D[i] < 0) psd = false;

    for (size_t j = i + 1; j < n; j++) {
      coeff sum2 = 0;
      for (size_t k = 0; k < i; k++)
        sum2 += L[j][k] * L[i][k] * D[k];
      coeff r = Q[j][i] - sum2;

      if (D[i] != 0) {
        L[j][i] = r / D[i];
        L[j][i].canonicalize();
      }
      else {
        if (r != 0) psd = false;
        L[j][i] = 0;
      }
    }
  }
  return psd;
}

bool
mvpolyT::exactify(bool dorefine) {
  if(!solution.solved)
    throw std::runtime_error("It is necessary to solve first");

  const unsigned int tol_digits = 100;
  mpz_class maxden = 1000000;

  if (dorefine) {
    refine(512, tol_digits, 60);
    if (solution.refined)
      mpz_ui_pow_ui(maxden.get_mpz_t(), 10, tol_digits / 2);
  }

  continued_fractions(maxden);
  projection();

  if (!verify_identity()) {
    if (verbose) std::cout << "[exactify] FAILED: f != z^T Q z\n";
    return false;
  }

  std::vector<std::vector<coeff>> L;
  std::vector<coeff> D;
  bool psd = LDLT_factor(gram_exact(), L, D);
  if (verbose)
    std::cout << "[exactify] identity f = z^T Q z holds exactly; Q is "
              << (psd ? "PSD  => exact certificate" : "NOT PSD  => no certificate") << '\n';
  return psd;
}

/*
 * Similar to Cholesky decomposition
 * Using the LDLT avoid square roots
 * and still in ℚ[X] 
 */
bool
mvpolyT::LDLT() {
  if(!solution.solved)
    throw std::runtime_error("It is necessary to solve first");
  if(solution.y_exact.empty())
    throw std::runtime_error("It is necessary to exactify first");

  std::vector<std::vector<coeff>> Q = gram_exact();
  std::vector<std::vector<coeff>> L;
  std::vector<coeff> D;
  const bool psd = LDLT_factor(Q, L, D);

  // dbg
  std::cout<<"Q =";
  for(unsigned int i = 0;i<mb.size; i++) {
    for(unsigned j = 0; j<mb.size; j++)
      std::cout<<"\t"<<Q[i][j];
    std::cout<<'\n';
  }

  std::cout<<"L=";
  for(unsigned int i = 0; i< mb.size;i++) {
    for(unsigned int j = 0; j< mb.size; j++)
      std::cout<<"\t"<<L[i][j];
    std::cout<<"\n";
  }
  std::cout<<"\n";

  std::cout<<"D=";
  for(unsigned int i=0;i<mb.size;i++)
    std::cout<<"\t"<<D[i];
  std::cout<<"\n";

  /*
   * TODO modularize
   * mvpolynomial sum of squares just for debuging
   * therefore print doubles
   */
  std::cout << "sos reconstruction of mvpolynomial\n";
  bool first = true;
  for (unsigned int i = 0; i < mb.size; i++) {
    double d_val = D[i].get_d();
    if (D[i] == 0)
      continue;

    if (!first)
      std::cout << " +\n";
    first = false;

    std::cout << "(" << d_val << ")*(";
    bool firstmon = true;
    for (unsigned int j = 0; j < mb.size; j++) {
      if (L[j][i] == 0)
        continue;

      if (!firstmon) {
        if (L[j][i] > 0)
          std::cout << " + ";
        else 
          std::cout << " - ";
      }
      bool firstcur = firstmon;
      firstmon = false;

      double l_val = L[j][i].get_d();

      if (firstcur) {
        if (l_val == -1.0)
          std::cout << "-";
        else if (l_val != 1.0)
          std::cout << l_val << "*";
      } else {
        if (std::abs(l_val) != 1.0)
          std::cout << std::abs(l_val) << "*";
      }

      bool vars = false;
      for (size_t v = 0; v < mb.monomials[j].size(); v++) {
        if (mb.monomials[j][v] == 0)
          continue;
        
        std::cout << "x_" << (v + 1);
        if (mb.monomials[j][v] != 1)
          std::cout << "^" << mb.monomials[j][v];
        std::cout << " ";
        vars = true;
      }
      if (!vars && (l_val == 1.0 || l_val == -1.0))
        std::cout << "1";
    }
    std::cout << ")^2";
  }
  std::cout << "\n";
  std::cout << (psd ? "Q is PSD (all D_i >= 0)\n" : "Q is NOT PSD\n");
  return psd;
}

/*
 *	Just for debuging
 */
void 
mvpolyT::dbg () {
  if (terms.empty()) {
    std::cout<<0<<'\n';
    return;
  }

  bool first_term = true;
  std::map<exps, coeff>::iterator it;
  for (it = terms.begin(); it != terms.end(); it++) {
    if (!first_term) {
      if(it->second > 0)
	std::cout<<" + ";
      else
	std::cout<<" - ";
    }
    else {
      /*
       * first_term
       */
      if(it->second < 0)
	std::cout<<" - ";
    }
    first_term = false;
   
    /*
     * gmpxx abs()
     */
    std::cout<<abs(it->second)<<" ";
    
    for (size_t i = 0; i<it->first.size(); i++) {
      if(it->first[i] == 0)
	continue;
      else
	std::cout<<" x_"<<i+1<<"^"<<it->first[i]<<" ";
    }
  }
  std::cout<<'\n';
}

void 
monomialbasisT::dbg() const {
  std::cout<<"Monomial basis\n------------------\n\n";
  std::cout<<"Cardinality: |m| = "<<size<<'\n';
  std::cout<<"Max-degree: d = "<<degree<<'\n'; 
  std::cout<<(homogeneous ? "\nThis is a base for a Form" : 
	      "\nThis is a base for a Polynomial")<<'\n';

  for (size_t i = 0; i<size; i++) {
    std::cout<<"m["<<i<<"] = [";
    for (size_t k = 0; k<monomials[i].size(); k++) {
      std::cout<<monomials[i][k];
      if (k < monomials[i].size() - 1)
	std::cout << ", ";
    }
    std::cout << "]";

    std::cout << " (";
    bool first = true;
    for (size_t k = 0; k < monomials[i].size(); k++) {
      if (monomials[i][k] > 0) {
	if(!first)
	  std::cout << " ";
	
	std::cout << "x_" << (k+1);

	if(monomials[i][k] !=1)
	  std::cout << "^" << monomials[i][k];

	first = false;
      } 
    }
    if (first)
      std::cout<<"1";
    std::cout<<")"<<'\n';
  }
}

void 
mvpolyT::mbdbg() {
  mb.dbg();
}
