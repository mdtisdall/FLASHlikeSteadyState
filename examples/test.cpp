#include "StructuredLUSolver.h"

#include "SequenceEventVectors.h"

#include "MPRAGE.h"

#include <vector>

#include <iostream>

typedef std::vector<float> vec_type;

typedef SequenceEventVectors<vec_type> SeqVecsT;

typedef MPRAGE<SeqVecsT> MPRAGET;

int main() {

  SeqVecsT seqVecs(1000);

  MPRAGET mprage(4, 400, 600, 176, 6); 
  
  mprage(&seqVecs);

  StructuredLUSolver<vec_type> solver(seqVecs.a.size());

  vec_type m(seqVecs.a.size());
  
  solver(&(seqVecs.a), &(seqVecs.b), &m);
  
  vec_type longitudinal(m.size());

  longitudinal[0] = m[m.size()-1];

  BLAS::copy(m.size() - 1, m.data(), 1, longitudinal.data() + 1, 1);

  vec_type trans;

  vec_type::iterator longIt = longitudinal.begin();

  mprage(&longIt, &trans);
  
  {
    std::cout << "trans: ";
    for(vec_type::iterator tIt = trans.begin(); tIt != trans.end(); tIt++) {
      std::cout << *tIt << ", ";
    }
    std::cout << std::endl;
  }
  
  return 0;
}
