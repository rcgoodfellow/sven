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

#ifdef DEBUG
  internal::RT::map_name(Q, "Q");
  internal::RT::map_name(H, "H");
  internal::RT::map_name(A, "A");
  internal::RT::map_name(x0, "x0");
  internal::RT::map_name(b, "b");
  internal::RT::map_name(r0, "r0");
#endif
}

size_t Arnoldi::N() const { return A.m(); }
size_t Arnoldi::n() const { return _n; }

using std::mutex;

void Arnoldi::operator()()
{
  //mutex print_mutex{};

  r0 = b - A*x0;
  r0 /= norm(r0);
  Q.C(0) = r0;

  //print_mutex.lock();
  //std::cout << "Q" << std::endl;
  //std::cout << Q << std::endl;
  //print_mutex.unlock();

  Vector x = Vector::Zero(N());
  Vector s = Vector::Zero(n());
  Scalar xn(0.0);
#ifdef DEBUG
  internal::RT::map_name(x, "x");
  internal::RT::map_name(s, "s");
  internal::RT::map_name(xn, "xn");
#endif

  for(size_t i=0; i<1; ++i)
  {
    //new subspace element
    Q.C(i+1) = A*Q.C(i);

    //Q.wait();
    MARKER("a"); 

    //ortho
    H.C(i) = T(Q.C(0,i)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i) * H.C(i);

    //H.wait();
    //Q.wait();
    MARKER("b");

    //reortho
    s = T(Q.C(0,i)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i) * s;
    H.C(i) += s;
    
    //normo
    x = Q.C(i+1);
    xn = norm(x);
    H.C(i)(i+1) = xn;
    Q.C(i+1) = x / xn;
    
    //H.wait();
    //Q.wait();
    MARKER("c");
  }

  /*
  while(internal::RT::get().active_jobs > 0 || !internal::RT::Q().empty())
  {
    usleep(10);
  }
    
  */

}

