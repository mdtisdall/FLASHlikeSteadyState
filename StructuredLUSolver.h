#ifndef StructuredLUSolver_h
#define StructuredLUSolver_h

#include "blas_local.h"

template <typename vec_type>
class StructuredLUSolver
{
  public:
    StructuredLUSolver (
      const size_t length
    ) :
      length(length)
     ,y(length)
     ,u(length)
    {}

    typedef typename vec_type::value_type value_type;

    void operator()(const vec_type *a, const vec_type *b, vec_type *m)
    {
      {
        // Solve Ly = b where L has 1 on the diagonal, and
        // a[1] ... a[length-1] on the subdiagonal
        y[0] = (*b)[0];
    
        u[0] = -(*a)[0];

        typename vec_type::iterator yIt = y.begin() + 1;
        typename vec_type::const_iterator yTrailingIt = y.begin();
        typename vec_type::iterator uIt = u.begin() + 1;
        typename vec_type::const_iterator uTrailingIt = u.begin();
        typename vec_type::const_iterator bIt = b->begin() + 1;
        typename vec_type::const_iterator aIt = a->begin() + 1;

        while(yIt != y.end()) {
          *yIt = *bIt - *yTrailingIt * *aIt;
          *uIt = *uTrailingIt * (-*aIt);
          yIt++;
          bIt++;
          yTrailingIt++;
          aIt++;
          uIt++;
          uTrailingIt++;
        }
      }

      {
        // Solve Um = y where U is the identity matrix with
        // u_i = prod_{i=0}^{length - 1} -a[i] subtracted from the
        // ith row of the last column. Since this is an atomic
        // upper triangular matrix, its inverse is constructed by negating
        // the off-diagonal elements.
        value_type lastM =
          (*m)[length-1] = y[length-1] / (value_type(1) - u[length-1]);

        // having computed the last entry in m, we can write every other
        // entry as the sum of two entries: the current y, and the current u
        // scaled by the last entry in m.
        BLAS::copy(length - 1, y.data(), 1, m->data(), 1);
        BLAS::axpy(length - 1, &lastM, u.data(), 1, m->data(), 1);
      }
    }

  protected:
    size_t length;
    vec_type y;
    vec_type u;
};

#endif
