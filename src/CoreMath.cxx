/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 15 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/CoreMath.hxx"

using namespace sven;
using std::max;
using std::initializer_list;
using std::move;
using std::unique_lock;
using std::lock_guard;
using std::mutex;
using std::condition_variable;
using std::runtime_error;
using std::string;
using std::to_string;

//~= Vector ~=-----------------------------------------------------------------
Vector::Vector(size_t n) 
  : _n{n}, _{alloc<double>(n)},
    _state{new ObjectState},
    _mtx{new mutex},
    _cnd{new condition_variable}
{ 
  *_state = ObjectState::Materializing;
}

Vector::Vector(initializer_list<double> xs)
  : Vector(xs.size())
{
  size_t i{0};
  for(auto d : xs){ _[i++] = d; }
  *_state = ObjectState::SolidState;
  _cnd->notify_all();
}

Vector Vector::Zero(size_t n)
{
  Vector x(n);
  //x._[0:n] = 0;
  for(size_t i=0; i<n; ++i)
  {
    x._[i] = 0;
  }
  *x._state = ObjectState::SolidState;
  x._cnd->notify_all();
  return x;
}

size_t Vector::n() const { return _n; }

double & Vector::operator()(size_t i)
{ 
  unique_lock<mutex> lk{*_mtx};
  switch(*_state)
  {
    case ObjectState::Materializing: 
      _cnd->wait(lk);
      lk.unlock();
      return _[i];

    case ObjectState::SolidState:    
      lk.unlock(); 
      return _[i];

    case ObjectState::Vapor:         
    default:                         
      lk.unlock();
      throw runtime_error("attempt to access vaporized object" + CRIME_SCENE);
  }
}

ObjectState Vector::state() const { return *_state; }

#include <iostream>

bool Vector::operator== (const Vector &x)
{
  unique_lock<mutex> lk_me{*_mtx}, lk_x{*x._mtx};
  if(*_state != ObjectState::SolidState) 
  { 
    lk_x.unlock();
    _cnd->wait(lk_me); 
    lk_x.lock();
  }
  if(*x._state != ObjectState::SolidState) 
  { 
    lk_me.unlock();
    x._cnd->wait(lk_x); 
    lk_me.lock();
  }

  if(x._n != _n){ return false; }

  for(size_t i=0; i<_n; ++i)
  {
    if(x._[i] != _[i]) 
    { 
      std::cout << i << std::endl;
      std::cout << x._[i] << std::endl;
      std::cout << _[i] << std::endl;
      return false; 
    }  
  }

  lk_me.unlock();
  lk_x.unlock();
  return true;
}

Vector sven::operator+ (const Vector &a, const Vector &b)
{
  if(a.n() != b.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector ab(a.n());
  internal::Thunk t = [a,b,ab](){ op_plus_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_plus_impl(const Vector a, const Vector b, Vector ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  for(size_t i=0; i<a.n(); ++i)
  { 
    //ab._[0:ab.n()] = a._[0:a.n()] + b._[0:b.n()]; 
    for(size_t i=0; i<ab.n(); ++i)
    {
      ab._[i] = a._[i] + b._[i]; 
    }
    *ab._state = ObjectState::SolidState;
    ab._cnd->notify_all();
  }
}

Vector sven::operator- (const Vector &a, const Vector &b)
{
  if(a.n() != b.n())
  {
    throw runtime_error("non-conformal operaton:" + CRIME_SCENE);
  }

  Vector ab(a.n());
  internal::Thunk t = [a,b,ab](){ op_sub_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_sub_impl(const Vector a, const Vector b, Vector ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  for(size_t i=0; i<a.n(); ++i)
  {
    //ab._[0:ab.n()] = a._[0:a.n()] - b._[0:b.n()]; 
    for(size_t i=0; i<ab.n(); ++i)
    {
      ab._[i] = a._[i] - b._[i]; 
    }
    *ab._state = ObjectState::SolidState;
    ab._cnd->notify_all();
  }
}


Scalar sven::operator* (const Vector &a, const Vector &b)
{
  if(a.n() != b.n()){ 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_dot_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

double sven::dot(size_t n, double *a, double *b)
{
  double dot{0};
  for(size_t i=0; i<n; ++i)
  {
    dot += a[i] * b[i];
  }
  return dot;
}

void sven::op_dot_impl(const Vector a, const Vector b, Scalar ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  *ab._ = dot(a.n(), a._, b._);
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

Vector sven::operator/ (const Vector &a, const Scalar &b)
{
  Vector ab(a.n());
  internal::Thunk t = [a,b,ab](){ op_div_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_div_impl(const Vector a, const Scalar b, Vector ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  //ab._[0:a._n] = a._[0:a._n] / *b._;
  for(size_t i=0; i<a._n; ++i)
  {
    ab._[i] = a._[i] / *b._;
  }
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

Vector sven::operator* (const Vector &a, const Scalar &b)
{
  Vector ab(a.n());
  internal::Thunk t = [a,b,ab](){ op_mul_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_mul_impl(const Vector a, const Scalar b, Vector ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  //ab._[0:a._n] = a._[0:a._n] * *b._;
  for(size_t i=0; i<a._n; ++i)
  {
    ab._[i] = a._[i] * *b._;
  }
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

//~= Scalar ~=-----------------------------------------------------------------
Scalar::Scalar() : _{alloc<double>(1)} 
{
  *_state = ObjectState::Materializing;
} 

Scalar::Scalar(double value) : Scalar() { *_ = value; }

double & Scalar::operator()() 
{ 
  unique_lock<mutex> lk{*_mtx};
  switch(*_state)
  {
    case ObjectState::Materializing:
      _cnd->wait(lk);
      lk.unlock();
      return *_;
      
    case ObjectState::SolidState:
      lk.unlock();
      return *_;

    case ObjectState::Vapor:
    default:
      lk.unlock();
      throw runtime_error("attempt to access vaporized object" + CRIME_SCENE);
  }
}

ObjectState Scalar::state() const { return *_state; }

Scalar sven::operator+ (const Scalar &a, const Scalar &b)
{
  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_plus_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_plus_impl(const Scalar a, const Scalar b, Scalar ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  *ab._ = *a._ + *b._;
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

Scalar sven::operator- (const Scalar &a, const Scalar &b)
{
  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_sub_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_sub_impl(const Scalar a, const Scalar b, Scalar ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  *ab._ = *a._ - *b._;
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

Scalar sven::operator* (const Scalar &a, const Scalar &b)
{
  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_mul_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_mul_impl(const Scalar a, const Scalar b, Scalar ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  *ab._ = *a._ * *b._;
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

Scalar sven::operator/ (const Scalar &a, const Scalar &b)
{
  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_div_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_div_impl(const Scalar a, const Scalar b, Scalar ab)
{
  lock_guard<mutex> lk_a{*a._mtx}, lk_b{*b._mtx}, lk_ab{*ab._mtx};
  *ab._ = *a._ / *b._;
  *ab._state = ObjectState::SolidState;
  ab._cnd->notify_all();
}

//~= Matrix ~=-----------------------------------------------------------------
Matrix::Matrix(size_t m, size_t n)
  : _m{m}, _n{n},
    _{alloc<double>(m*n)}
{
  *_state = ObjectState::Materializing;
  _cnd->notify_all();
}

Matrix Matrix::Zero(size_t m, size_t n)
{
  Matrix A{m,n};
  //A._[0:m*n] = 0;
  for(size_t i=0; i<m*n; ++i)
  {
    A._[i] = 0;
  }
  *A._state = ObjectState::SolidState;
  A._cnd->notify_all();
  return A;
}

Matrix Matrix::Identity(size_t m, size_t n)
{
  Matrix A = Matrix::Zero(m,n);
  for(size_t i=0; i<max(m,n); ++i){ A._[i*n + i] = 1; }
  *A._state = ObjectState::SolidState;
  A._cnd->notify_all();
  return A;
}

size_t Matrix::m() const { return _m; }
size_t Matrix::n() const { return _n; }

Vector sven::operator* (const Matrix &A, const Vector &x)
{
  if(A.n() != x.n()){ 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector Ax(x.n());
  internal::Thunk t = [A,x,Ax](){ op_mul_impl(A,x,Ax); };
  internal::RT::Q().push(t);
  return Ax;
}

void sven::multi_dot(size_t n, double *a, double *b, double *ab, 
    CountdownLatch &cl)
{
  *ab = dot(n, a, b);
  --cl;
}

void sven::op_mul_impl(const Matrix A, const Vector x, Vector Ax)
{
  lock_guard<mutex> lk_A{*A._mtx}, lk_x{*x._mtx}, lk_Ax{*Ax._mtx};
  CountdownLatch cl{static_cast<int>(A.n())};
  for(size_t i=0; i<A.n(); ++i)
  {
    //Ax._[i] = dot(A.n(), &A._[i*A.n()], x._);
    internal::Thunk t = [A, x, Ax, &cl, i]()
    { 
      multi_dot(A.n(), &A._[i*A.n()], x._, &Ax._[i], cl); 
    };

    internal::RT::Q().push(t);
  }
  cl.wait();
  *Ax._state = ObjectState::SolidState;
  Ax._cnd->notify_all();
}

Matrix sven::operator* (const Matrix &A, const Matrix &B)
{
  if(A.n() != B.m()){ 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Matrix AB(A.m(), B.n());
  internal::Thunk t = [A,B,AB](){ op_mul_impl(A,B,AB); };
  internal::RT::Q().push(t);
  return AB;

}

void sven::op_mul_impl(const Matrix A, const Matrix B, Matrix AB)
{
  //TODO: You are here!
}
