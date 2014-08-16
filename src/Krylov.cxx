/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 16 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/Krylov.hxx"

using namespace sven;

Arnoldi::Arnoldi(size_t N, size_t n, Vector x0)
  : Q{N,n}, H{n,n},
    x0{x0},
    _N{N}, _n{n}
{}

size_t Arnoldi::N() const { return _N; }
size_t Arnoldi::n() const { return _n; }
