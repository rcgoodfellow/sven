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
    Arnoldi(size_t N, size_t n, Vector x0);

    Matrix Q, H;
    Vector x0;

    size_t N() const, n() const;

  private:
    size_t _N, _n;

};

}

#endif
