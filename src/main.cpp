/*
 *	Converts an multivariate polynomial sparse point value representation to the sparse sdpa format
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


  // TODO projection in affine subspace
  F.exactify();

  F.LDLT();

  // just for debuging
  F.mbdbg();
  }
  catch (const std::exception& e){
    std::cerr << "Error: "<<e.what()<<'\n';
    return(1);
  }

  return(0);
}
