/*
 *	Reads a multivariate polynomial (sparse representation) from stdin,
 *	solves the Gram-matrix SDP with SCIP-SDP, refines + exactifies the
 *	solution and verifies the exact SOS certificate.
 */
#include <exception>
#include <iostream>
#include "../include/mvpolynomial.hpp"

int main () {
  
  try{
  /*
   * Reading the multivariate polynomial here
   */
  mvpolyT F; //(x_1, ..., x_n)
  
  std::cout<<"Polynomial"<<'\n';
  F.dbg();

  std::cout << '\n';

  /*
   * Solve using the SCIP-SDP API
   * and the SCIPpp wrapper.
   */
  F.solvesdp();

  /*
   * refine (Gauss-Newton) + continued fractions + projection
   * + exact verification (f == z^T Q z and Q PSD)
   */
  bool certified = F.exactify();

  F.LDLT();

  // just for debuging
  F.mbdbg();

  if (!certified) {
    std::cerr << "NOT certified: f is not strictly inside the SOS cone "
                 "(or the SDP solution is too inaccurate).\n";
    return 2;
  }
  std::cout << "CERTIFIED: exact SOS decomposition in Q.\n";
  }
  catch (const std::exception& e){
    std::cerr << "Error: "<<e.what()<<'\n';
    return(1);
  }

  return(0);
}
