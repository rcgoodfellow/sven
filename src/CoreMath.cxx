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
using std::vector;

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

  bool result{true};
  for(size_t i=0; i<_n; ++i)
  {
    if(x._[i] != _[i]) { result = false; break; }  
  }

  lk_me.unlock();
  lk_x.unlock();
  return result;
}
    
std::ostream & sven::operator<< (std::ostream &o, Vector &x)
{
  for(size_t i=0; i<x._n; ++i)
  {
    o << x._[i] << " ";
  }
  return o;
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

double sven::dot(size_t n, 
    double *a, size_t a_stride, 
    double *b, size_t b_stride)
{
  double d{0};
  for(size_t i=0; i<n; ++i)
  {
    d += a[i*a_stride] * b[i*b_stride];
  }
  return d;
}

double sven::dot(size_t n, double *a, double *b)
{
  double d{0};
  for(size_t i=0; i<n; ++i)
  {
    d += a[i] * b[i];
  }
  return d;
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
}

Matrix::Matrix(size_t m, size_t n, vector<double> values)
  : Matrix(m, n)
{
  if(values.size() != m*n)
  {
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  for(size_t i=0; i<m*n; ++i) { _[i] = values[i]; }
  *_state = ObjectState::SolidState;
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

bool Matrix::operator== (const Matrix &A)
{
  unique_lock<mutex> lk_me{*_mtx}, lk_A{*A._mtx};
  if(*_state != ObjectState::SolidState)
  {
    lk_A.unlock();
    _cnd->wait(lk_me);
    lk_A.lock();
  }
  if(*A._state != ObjectState::SolidState)
  {
    lk_me.unlock();
    A._cnd->wait(lk_A);
    lk_me.lock();
  }

  if(A._m != _m || A._n != _n){ return false; }

  bool result{true};
  for(size_t i=0; i<_n*_m; ++i)
  {
    if(A._[i] != _[i]){ result = false; break; }
  }

  lk_me.unlock();
  lk_A.unlock();
  return result;
}

double & Matrix::operator()(size_t i, size_t j)
{
  unique_lock<mutex> lk{*_mtx};
  switch(*_state)
  {
    case ObjectState::Materializing: 
      _cnd->wait(lk);
      lk.unlock();
      return _[i*_n + j];

    case ObjectState::SolidState:    
      lk.unlock(); 
      return _[i*_n + j];

    case ObjectState::Vapor:         
    default:                         
      lk.unlock();
      throw runtime_error("attempt to access vaporized object" + CRIME_SCENE);
  }
}

size_t Matrix::m() const { return _m; }
size_t Matrix::n() const { return _n; }

Vector sven::operator* (const Matrix &A, const Vector &x)
{
  if(A.n() != x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector Ax(x.n());
  internal::Thunk t = [A,x,Ax](){ op_mul_impl(A,x,Ax); };
  internal::RT::Q().push(t);
  return Ax;
}

void sven::multi_dot(size_t n, 
    double *a, size_t a_stride,
    double *b, size_t b_stride,
    double *ab, 
    CountdownLatch &cl)
{
  *ab = dot(n, a, a_stride, b, b_stride);
  --cl;
}

void sven::multi_dot(size_t n, 
    double *a,
    double *b,
    double *ab, 
    CountdownLatch &cl)
{
  multi_dot(n, a, 1, b, 1, ab, cl);
}

void sven::op_mul_impl(const Matrix A, const Vector x, Vector Ax)
{
  internal::Thunk th = [A,x,Ax]()
  {
    lock_guard<mutex> lk_A{*A._mtx}, lk_x{*x._mtx}, lk_Ax{*Ax._mtx};

    CountdownLatch cl{static_cast<int>(A.n())};
    for(size_t i=0; i<A.n(); ++i)
    {
      internal::Thunk t = [A, x, Ax, &cl, i]()
      { 
        multi_dot(A.n(), &A._[i*A.n()], x._, &Ax._[i], cl); 
      };

      internal::RT::Q().push(t);
    }
    cl.wait();
    *Ax._state = ObjectState::SolidState;
    Ax._cnd->notify_all();
  };

  internal::RT::Q().push(th);
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
  internal::Thunk th = [A,B,AB]()
  {
    lock_guard<mutex> lk_A{*A._mtx}, lk_B{*B._mtx}, lk_AB{*AB._mtx};
    CountdownLatch cl{static_cast<int>(A._m*B._n)};
    for(size_t i=0; i<A._m; ++i)
    {
      for(size_t j=0; j<B._n; ++j)
      {
        internal::Thunk t = [A, B, AB, &cl, i, j]()
        {
          multi_dot(A._n, 
              &A._[i*A._n], 1, 
              &B._[j], B._n, 
              &AB._[i*A._n + j], cl);
        };

        internal::RT::Q().push(t);
      }
    }
    cl.wait();
    *AB._state = ObjectState::SolidState;
    AB._cnd->notify_all();
  };

  internal::RT::Q().push(th);
}

Column Matrix::C(size_t index)
{
  return Column(this, index);
}

//~=~ Column ~=~---------------------------------------------------------------

Column::Column(Matrix *origin, size_t index)
  : _origin{origin}, _index{index}
{}
    
Column & Column::operator= (const Vector &x)
{

  unique_lock<mutex> lk{*x._mtx};
  if(_origin->_m != x._n)
  { 
    lk.unlock();
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }
  if(*x._state != ObjectState::SolidState) { x._cnd->wait(lk); }

  for(size_t i=0; i<_origin->_m; ++i)
  {
    _origin->_[i*_origin->_n + _index] = x._[i]; 
  }
  lk.unlock();
  return *this;
}

bool Column::operator== (const Vector &x)
{
  unique_lock<mutex> lk{*x._mtx};
  if(*x._state != ObjectState::SolidState) { x._cnd->wait(lk); }
  
  if(_origin->_m != x._n) { lk.unlock(); return false; }

  bool result{true};
  for(size_t i=0; i<_origin->_m; ++i)
  {
    if(_origin->_[i*_origin->_n + _index] != x._[i]){ result = false; break; }
  }
  lk.unlock();
  return result;
}

//~=~ SparseMatrix ~=~---------------------------------------------------------

SparseMatrix::SparseMatrix(size_t m, size_t n, size_t z)
  : _m{m}, _n{n}, _z{z},
    _r{alloc<size_t>(m)},
    _c{alloc<size_t>(m*z)},
    _v{alloc<double>(m*z)}
{
  *_state = ObjectState::Materializing;
}
    
SparseMatrix::SparseMatrix(size_t m, size_t n, size_t z, 
    std::vector<size_t> rs, std::vector<size_t> cs, 
    std::vector<double> vs)
  : SparseMatrix(m, n, z)
{
  if(rs.size() > m)
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }
  if(cs.size() > m*n)
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }
  if(vs.size() > m*n)
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }
  if(vs.size() != cs.size())
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }

  size_t j=0;
  for(size_t i=0; i<rs.size(); ++i)
  {
    if(rs[i] > z)
    {
      throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
    }
    _r[i] = rs[i];
    for(size_t k=0; k<rs[i]; ++k, ++j)
    {
      if(cs[j] > n)
      {
        throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
      }
      _c[i*z + k] = cs[j];
      _v[i*z + k] = vs[j];
    }
  }

  if(j != vs.size())
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }
}

SparseMatrix SparseMatrix::Identity(size_t m, size_t n, size_t z)
{
  SparseMatrix A(m, n, z);
  for(size_t i=0; i<max(m,n); ++i)
  {
    A._r[i] = 1;
    A._c[i*z] = i;
    A._v[i*z] = 1;
  }
  *A._state = ObjectState::SolidState;
  A._cnd->notify_all();
  return A;
}

double & SparseMatrix::_at(size_t i, size_t j)
{
  for(size_t k=0; k<_r[i]; ++k)
  {
    if(_c[_z*i + k] == j){return _v[_z*i + k];}
  }

  throw runtime_error("unreal sparse matrix coordinates" + CRIME_SCENE);
}

double & SparseMatrix::operator() (size_t i, size_t j)
{
  if(i >= _m || j >= _n)
  { 
    throw runtime_error("unreal sparse matrix coordinates" + CRIME_SCENE);
  }
  unique_lock<mutex> lk{*_mtx};
  switch(*_state)
  {
    case ObjectState::Materializing: 
      _cnd->wait(lk);
      //fallthrough

    case ObjectState::SolidState:    
      lk.unlock(); 
      return _at(i, j);

    case ObjectState::Vapor:         
    default:                         
      lk.unlock();
      throw runtime_error("attempt to access vaporized object" + CRIME_SCENE);
  }
  
}

size_t SparseMatrix::m() const { return _m; }
size_t SparseMatrix::n() const { return _n; }
size_t SparseMatrix::z() const { return _z; }

Vector sven::operator* (const SparseMatrix &A, const Vector &x)
{
  if(A.n() != x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector Ax(x.n());
  internal::Thunk t = [A,x,Ax](){ op_mul_impl(A,x,Ax); };
  internal::RT::Q().push(t);
  return Ax;
}

void sven::multi_sparse_dot(size_t rz,
    size_t *A_c, double *A_v,
    double *x,
    double *Ax,
    CountdownLatch &cl)
{
  double d{0};
  for(size_t i=0; i<rz; ++i)
  {
    d += A_v[i] * x[A_c[i]];
  }
  *Ax = d;
  --cl;
}

void sven::op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax)
{
  internal::Thunk th = [A,x,Ax]()
  {
    lock_guard<mutex> lk_A{*A._mtx}, lk_x{*x._mtx}, lk_Ax{*Ax._mtx};

    CountdownLatch cl{static_cast<int>(A.n())};
    for(size_t i=0; i<A.n(); ++i)
    {
      internal::Thunk t = [A,x,Ax,&cl,i]()
      {
        multi_sparse_dot(A._r[i], &A._c[i*A._z], &A._v[i*A._z], x._, &Ax._[i], 
            cl);
      };
      internal::RT::Q().push(t);
    }
    cl.wait();
    *Ax._state = ObjectState::SolidState;
    Ax._cnd->notify_all();
  };
  
  internal::RT::Q().push(th);
}
