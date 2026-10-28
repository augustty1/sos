#ifndef MVPOLYNOMIAL_HPP
#define MVPOLYNOMIAL_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <gmpxx.h>

using exps = std::vector<unsigned int>;

/*Polynomials are in ℚ[X]*/
using coeff = mpq_class;

struct constraintT {
  std::map<unsigned int, coeff> lhs;
  coeff rhs;
};

/*
* SCIP-SDP solution is in floating point representation
* eventually, the continued fraction will transform
* it from double to coeffs = mpq_class again.
*
* Pipeline: y(double) --refine--> y_r (mpf, ~150 digits)
*           --continued fractions--> y_c (mpq) --projection--> y_exact (in the affine subspace)
* 
* All three vectors use the same indexing. k runs over the upper triangle
* (i <= j) of the Gram matrix, row by row.
*
*/
struct sdpsolT {
  std::vector<double> y;
  std::vector<mpf_class> y_refined;
  std::vector<coeff> y_exact;
  unsigned int nvars;
  size_t size;
  bool refined;
  bool solved;

  sdpsolT() : nvars(0), size(0),refined(false), solved(false) {}
};

// TODO Newton Polytope optimization
class monomialbasisT {
  public:
  std::vector<exps> monomials;
  size_t size;
  bool homogeneous;
  
  /*half of max degree of mvpolyT*/
  unsigned int degree;

  void dbg() const;
};

class mvpolyT {
  unsigned int nvars;
  std::map<exps, coeff> terms;
  monomialbasisT mb;
  bool ismbcomputed;
  sdpsolT solution;  

  unsigned int degree() const;
  bool homogeneous() const;
  void rec_monomialbasis(unsigned int v, unsigned int remain, exps& curr, std::vector<exps>& ans) const;

  /*
   * Exatification methods
   */
  void continued_fractions(const mpz_class& maxden);
  void refine(unsigned int prec_bits = 512, unsigned int tol_digits = 100, unsigned int maxiter=60);
  void projection();

  /*
   * Certificate checks
   * all in exact arithmetic
   */
  bool verify_identity() const;
  bool LDLT_factor(const std::vector<std::vector<coeff>>& Q,
                    std::vector<std::vector<coeff>>& L,
                    std::vector<coeff>& D) const;
  std::vector<std::vector<coeff>> gram_exact() const;

  public:
  bool verbose;

  mvpolyT ();
  void addterm(const exps& _e, coeff _c);
  const monomialbasisT& monomialbasis();
  std::vector<constraintT> gram();

/*
 * In this branch, this method uses the SCIP-SDP API
 * and the SCIPpp wrapper.
 */
  void solvesdp();

  /*
   * Call the private exatification methods
   */
  bool exactify(bool dorefine = true);  
  
  /*
   * matrix decomposition
   * returns true iff Q is PSD
   */
  bool LDLT();
 
  void dbg();
  void mbdbg();
};

#endif // MVPOLYNOMIAL_HPP
