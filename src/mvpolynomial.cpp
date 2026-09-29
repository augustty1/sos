#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "../include/mvpolynomial.hpp"

mvpolyT::mvpolyT () {
  std::cin>>nvars;
  double temp_c;
  exps e(nvars);
  ismbcomputed = false;
  
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
 * constraints. This constraints will be translated to the .dat-s format
 * for a Semidefinite programming solver.
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
	coeff c = (i==j)? 1.0 : 2.0;

	monomialindex[product].push_back({y_k, c});
      }
    }

    // Generating linear constraints (new)
    // Check all possible monimials in momnomilindex and compair with the polynomial
    for (auto const& pair : monomialindex) {
        constraintT ctr;
        auto term_it = terms.find(pair.first);
    ctr.rhs = (term_it != terms.end()) ? term_it->second : coeff(0);

    const std::vector<std::pair<unsigned int, coeff>>& indexandcoeffs = pair.second;
    
    for(size_t v = 0; v < indexandcoeffs.size(); v++) {
        unsigned int index = indexandcoeffs[v].first;
        coeff c = indexandcoeffs[v].second;
        ctr.lhs[index] = c;
    }
    
    ans.push_back(ctr);

    }
    
    return(ans);
}

/*
 * .dat-s is the general format for sparse mixed semidefinite programming
 * [more info](https://github.com/scipopt/SCIP-SDP/blob/main/sdpa_format.txt)
 */
void 
mvpolyT::todats(const std::string& filename) {
  const monomialbasisT& mb = monomialbasis();
  std::vector<constraintT>ctrs = gram();
  
  unsigned int numvars = (mb.size*(mb.size + 1)) / 2;
  /*
   * To gen == contraint, we need the >= and the <=
   */
  unsigned int numlps = 2*ctrs.size();

  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Couldn't open: "+ filename);

  file<<numvars<<'\n';
  file<<"2\n"; /*number of blocks sdp and lp constraints*/
  file<<mb.size<<" -"<<numlps<<'\n'; /*size of the blocks - is for lp constraints*/
  for(unsigned int i=0; i<numvars;i++) /*fesiability*/
    file<<(i<numvars-1?"0 ":"0");
  file<<'\n';

  /*
   * sdp block
   * the sdp matrix is 1 based
   */
  unsigned int y_k = 1;
  for (unsigned int i = 0; i<mb.size; i++) {
    for(unsigned int j=i; j<mb.size; j++) {
      file<<y_k<<" 1 "<<(i+1)<<" "<<(j+1)<<" 1.0\n";
      y_k++;
    }
  }

  /*
   * lp constraints block
   * the lp rows is 1 based
   */
  unsigned int row=1;
  for(size_t ctr=0; ctr< ctrs.size(); ctr++) {
   file <<"0 2 "<<row<<" "<<row<<" "<<ctrs[ctr].rhs<<'\n';/*constant term*/

   std::map<unsigned int, coeff>::iterator it;
   for(it = ctrs[ctr].lhs.begin() ; it != ctrs[ctr].lhs.end(); it++) {
     file<<it->first<<" 2 "<<row<<" "<<row<<" "<<it->second<<'\n'; 
   }
   row++;

   /*sign inversion*/
    file<<"0 2 "<<row<<" "<<row<<" "<<-(ctrs[ctr].rhs)<<'\n';
    for(it = ctrs[ctr].lhs.begin(); it != ctrs[ctr].lhs.end(); it++) {
      file<<it->first<<" 2 "<<row<<" "<<row<<" "<< -(it->second)<<'\n';
    }
    row++;
  }
  file.close();
}

void 
mvpolyT::scipsdp(const std::string& filename) const {
  std::string command = "scipsdp <" + filename;
  int res = std::system(command.c_str());
  if(res !=0)
    std::cerr<<"\n Error in scipsdp"<<res<<'\n';
}

/*
 * After SCIP Solver we need to read data and transform again into a polynomial
 */
void 
mvpolyT::readsolution(const std::string& filename) {
  if(!ismbcomputed)
    throw std::runtime_error("It is necessary to compute monomial basis first");

  std::ifstream file(filename);
  std::string line;
  bool foundopt = false;
  unsigned int expect = (mb.size *(mb.size + 1))/2;
  solution.y.assign(expect,0.0);

  while(std::getline(file,line)) {
    if(line.find("optimal solution found") != std::string::npos)
      foundopt = true;

    if(line.find("x_") != std::string::npos && line.find("(obj:") != std::string::npos) {
      std::stringstream ss(line);
      std::string x;
      /*
       * SCIP -> double -> mpq_class -> LDLT()
       */
      //double temp_c;
      //ss>>x>>temp_c;

      //coeff c(temp_c);
      //c.canonicalize();
      double c;
      ss >> x >> c;

      size_t pos = x.find("_");
      if (pos!=std::string::npos) {
	        unsigned int i = std::atoi(x.c_str() + pos + 1);
	 if(i < expect)
	    solution.y[i] = c;
      }
    }
  }

  file.close();

  // another debug
  if(!foundopt)
    throw std::runtime_error("Optimal solution was not found");

  solution.nvars = solution.y.size();
  solution.size = mb.size;
  solution.solved = true;
/*
  for(int i = 0; i<solution.y.size(); i++)
    std::cout<<solution.y[i]<<'\n';
*/
}

/*
* Approximate a double with finite 
* continued fraction i.e. a better rational
* approximation: a/b with b > 0 there is no 
* rational number with smaller denominator 
* which is closer to the "real" number
*/
void
mvpolyT::continued_fractions(long maxden) {
    solution.y_exact.resize(solution.y.size());

    for (size_t k = 0; k < solution.y.size(); k++) {
        double x = solution.y[k];
        long sign = (x < 0) ? -1 : 1;
        x = std::abs(x);

        long long h0 = 0, h1 = 1;
        long long g0 = 1, g1 = 0;
        double y = x;
        long long a = static_cast<long long>(std::floor(y));

        long long h = a*h1 + h0;
        long long g = a*g1 + g0;

        while (g <= maxden && y!= static_cast<double>(a)) {
            double diff = y - static_cast<double>(a);
            if (diff == 0.0) break;
            y = 1.0 / diff;
            h0 = h1;
            h1 = h;
            g0 = g1;
            g1 = g;
            a = static_cast<long long>(std::floor(y));
            h = a*h1 + h0;
            g = a*g1 + g0;
        }
        if ( g > maxden ) {
            h = h1;
            g = g1;
        }

        mpq_class rational(static_cast<long>(sign*h), static_cast<unsigned long>(g));
        rational.canonicalize();
        solution.y_exact[k] = rational;

    }
    /*
     * Evaluate tau is compute the euclidian distance: 
     * tau^2 = sum_(i,j)(Q_i,j - ~Q_i,j)^2
     */
}

/*
 * Projection into the affine subspace of 
 * valid gram matrices.
 */
// TODO modulo check solution
// TODO data structure for monomial hash table
void
mvpolyT::projection() {
 if(solution.y_exact.empty()) {
       throw std::runtime_error("It is necessary to solve first");
   }
   
    std::cout<<"Projection\n------------\n\n"; /*dbg*/

    /*
     * expoents |-> y index and coefficient
     */
    std::map<exps, std::vector<std::pair<unsigned int, coeff>>> monomialindex;
    unsigned int k = 0;
    for (unsigned int i = 0; i < mb.size; i++) {
        for (unsigned int j = i; j < mb.size; j++) {
            exps product(nvars);
            for (unsigned int v = 0; v < nvars; v++) {
                product[v] = mb.monomials[i][v] + mb.monomials[j][v];
            }
            coeff c = (i == j) ? 1 : 2;
            monomialindex[product].push_back({k, c});
            k++;
        }
    }

    for (const auto& item : monomialindex) {
        /*
         * gamma = alpha + beta (tuples of expoents) 
         * entries that contribute to the monomial
         * in final polynomial
         */
        const exps& gamma = item.first;
        const auto& entries = item.second;

        /*
         * projection divisor
         */
        coeff n_gamma = 0;
        coeff current_val = 0;

        for (const auto& entry : entries) {
            unsigned int k = entry.first;
            coeff mult = entry.second;

            n_gamma += mult; 
            current_val += mult * solution.y_exact[k]; 
        }

        /*
         * At original polynomial
         */
        coeff target_val = 0;
        auto target_it = terms.find(gamma);
        if (target_it != terms.end()) {
            target_val = target_it->second;
        }
        coeff error = current_val - target_val;
        
        /*
         *  dbg
         */
        std::cout<<error<<'\n';

        /*
         * Apply correction
         */
        if (error != 0) {
            coeff correction = error / n_gamma;
            for (const auto& entry : entries) {
                unsigned int k = entry.first;
                solution.y_exact[k] -= correction;
                solution.y_exact[k].canonicalize(); 
            }
        }
    }
}

void
mvpolyT::exactify() {
   if(!solution.solved) {
       throw std::runtime_error("It is necessary to solve first");
   }

   continued_fractions(1e6);
   projection();
}

/*
 * Similar to Cholesky decomposition
 * Using the LDLT avoid square roots
 * and still in ℚ[X] 
 */
void 
mvpolyT::LDLT() {
  if(!solution.solved)
    throw std::runtime_error("It is necessary to solve first");

  /*
   * Gram matrix reconstruction
   * from solution
   */
  std::vector<std::vector<coeff>> Q( mb.size , std::vector<coeff>(mb.size, 0));
  unsigned k=0;
  for (unsigned int i = 0; i<mb.size ; i++){
    for (unsigned int j = i; j<mb.size; j++){
      // coeff q = solution.y[k];
      coeff q = solution.y_exact[k];
      // coeff q(solution.y[k]);
      q.canonicalize();
      Q[i][j] = q;
      Q[j][i] = q;
      k++;
    }
  }
    // dbg
    std::cout<<"Q =";
    for(unsigned int i = 0;i<mb.size; i++) {
      for(unsigned j = 0; j<mb.size; j++) {
	std::cout<<"\t"<<Q[i][j];
      }
      std::cout<<'\n';
    }
  /*
   * LDLT
   */
  std::vector<std::vector<coeff>> L(mb.size, std::vector<coeff>(mb.size, 0));
  std::vector<coeff> D(mb.size, 0);

  for (unsigned int i = 0; i<mb.size; i++)
    L[i][i]=1;
  for (unsigned int i = 0; i<mb.size; i++) {
    coeff sum = 0;
    for (unsigned int j=0; j<i; j++)
      sum += L[i][j]*L[i][j]*D[j];
    D[i] = Q[i][i] - sum;
    
    for (unsigned int j=i+1;j<mb.size;j++){
      coeff sum2 = 0;
      for(unsigned int k = 0; k<i; k++)
	sum2 += L[j][k] * L[i][k] * D[k];

      if (D[i] != 0) {
	    L[j][i] = (Q[j][i] - sum2) / D[i];
        L[j][i].canonicalize();
      }
      else {
       L[j][i]=0;
      }
    }
  }
  
  // TODO modularize
  // dbg
  
  std::cout<<"L=";
  for(unsigned int i = 0; i< mb.size;i++) {
    for(unsigned int j = 0; j< mb.size; j++) {
      std::cout<<"\t"<<L[i][j];
    }
    std::cout<<"\n";
  }
  std::cout<<"\n";

  std::cout<<"D=";
  for(unsigned int i=0;i<mb.size;i++) {
    std::cout<<"\t"<<D[i];
  }
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
    if (std::abs(d_val) <= 1e-10)
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
        if (mb.monomials[j][v] != 1) {
          std::cout << "^" << mb.monomials[j][v];
        }
        std::cout << " ";
        vars = true;
      }
      if (!vars && (l_val == 1.0 || l_val == -1.0)) {
        std::cout << "1";
      }
    }
    std::cout << ")^2";
  }
  std::cout << "\n";

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
