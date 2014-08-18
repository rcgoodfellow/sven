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

using std::mutex;

void Arnoldi::operator()()
{
  mutex print_mutex{};

  r0 = b - A*x0;
  r0 /= norm(r0);
  Q.C(0) = r0;

  print_mutex.lock();
  std::cout << "Q" << std::endl;
  std::cout << Q << std::endl;
  print_mutex.unlock();


  for(size_t i=0; i<_n-1; ++i)
  {
    Q.C(i+1) = A*Q.C(i);

    print_mutex.lock();
    std::cout << "Q1" << std::endl;
    std::cout << Q << std::endl;
    print_mutex.unlock();

    H.C(i) = T(Q.C(0,i)) * Q.C(i+1);
    std::cout << "H" << std::endl;
    std::cout << H << std::endl;

    Q.C(i+1) -= Q.C(0,i) * H.C(i);

    print_mutex.lock();
    std::cout << "Q2" << std::endl;
    std::cout << Q << std::endl;
    print_mutex.unlock();

    //reortho
    Vector s = T(Q.C(0,i)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i) * s;
    H.C(i) += s;
    Vector x = Q.C(i+1);
    Scalar xn = norm(x);
    xn.wait();
    H.C(i)(i+1) = xn();
    Q.C(i+1) = x / xn;

    
    print_mutex.lock();
    std::cout << "Q3" << std::endl;
    std::cout << Q << std::endl;
    print_mutex.unlock();
    
    print_mutex.lock();
    std::cout << "H" << std::endl;
    std::cout << H << std::endl;
    print_mutex.unlock();
  }

  while(internal::RT::get().active_jobs > 0 || !internal::RT::Q().empty())
  {
    usleep(10);
  }

}

