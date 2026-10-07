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
*/
struct sdpsolT {
  // std::vector<coeff> y;
  std::vector<double> y;
  std::vector<coeff> y_exact;
  unsigned int nvars;
  size_t size;
  bool solved;

  sdpsolT() : nvars(0), size(0), solved(false) {}
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
  void rec_monomialbasis(unsigned int v, unsigned int remain, 
       exps& curr, std::vector<exps>& ans) const;

  /*
   * Exatification methods
   */
  void continued_fractions(long maxden);
  void projection();
  

  public:
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
  void exactify();  
  
  /*
   * matrix decomposition
   */
  void LDLT();
 
  void dbg();
  void mbdbg();
};

#endif // MVPOLYNOMIAL_HPP
