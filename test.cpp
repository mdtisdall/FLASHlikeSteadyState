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
  
  {
    std::cout << "a: ";
    for(
      typename vec_type::iterator aIt =
        seqVecs.a.begin(); aIt != seqVecs.a.end(); aIt++)
    {
      std::cout << *aIt << ", ";
    }
    std::cout << std::endl;
  }
  
  {
    std::cout << "b: ";
    for(
      typename vec_type::iterator bIt =
        seqVecs.b.begin(); bIt != seqVecs.b.end(); bIt++)
    {
      std::cout << *bIt << ", ";
    }
    std::cout << std::endl;
  }

  solver(&(seqVecs.a), &(seqVecs.b), &m);

  std::cout << "m: ";
  for(vec_type::iterator mIt = m.begin(); mIt != m.end(); mIt++) {
    std::cout << *mIt << ", ";
  }
  std::cout << std::endl;
  
  return 0;
}
