/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 16 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/Krylov.hxx"

using namespace sven;
using std::runtime_error;

Arnoldi::Arnoldi(size_t n, SparseMatrix A, Vector x0, Vector b)
  : Q{Matrix::Zero(A.m(),n)}, 
    H{Matrix::Zero(n,n)}, 
    A{A}, x0{x0}, b{b}, r0(A.m()),
    _n{n}
{
  if(A.n() != x0.n() || A.n() != b.n())
  {
    throw runtime_error("non-conformal operation" + CRIME_SCENE);
  }
}

size_t Arnoldi::N() const { return A.m(); }
size_t Arnoldi::n() const { return _n; }

void Arnoldi::operator()()
{
  r0 = b - A*x0;
  r0 /= norm(r0);
  Q.C(0) = r0;

  for(size_t i=0; i<_n; ++i)
  {
    Q.C(i+1) = A*Q.C(i);
    std::cout << "Q" << std::endl;
    std::cout << Q << std::endl;

    H.C(i) = T(Q.C(0,i)) * Q.C(i+1);
    std::cout << "H" << std::endl;
    std::cout << H << std::endl;
  }

}

