#ifndef SVEN_KRYLOV_HXX
#define SVEN_KRYLOV_HXX

/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * The Krylov header contains template, class and function declarations for
 * doing Krylov subspace type things
 *
 * 16 August 2014
 * ~ ry
 *
 *****************************************************************************/

#include "CoreMath.hxx"

namespace sven
{

class Arnoldi
{
  public:
    Arnoldi() = delete;
    Arnoldi(size_t n, SparseMatrix A, Vector x0, Vector b);

    Matrix Q, H;
    SparseMatrix A;
    Vector x0, b;

    size_t N() const, n() const;
    void operator()();

    Vector r0;

  private:
    size_t _n;

};

}

#endif
