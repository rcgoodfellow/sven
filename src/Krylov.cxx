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
  : Q(A.m(),n), H(n,n), A{A}, x0{x0}, b{b}, r0(A.m()),
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
}

