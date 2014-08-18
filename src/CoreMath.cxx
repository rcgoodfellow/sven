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
using std::setprecision;
using std::fixed;
using std::adopt_lock;

//~= Vector ~=-----------------------------------------------------------------
Vector::Vector(size_t n) 
  : Object{alloc<double>(n)},
    _n{n}
{ 
  tick();
}

Vector::Vector(initializer_list<double> xs)
  : Vector(xs.size())
{
  size_t i{0};
  for(auto d : xs){ _data[i++] = d; }
  tock();
}

Vector::Vector(const Vector &x)
  : Object(x), _n{x.n()}
{ }

Vector::Vector(const Column &c)
  : Vector(c.origin->m())
{
  lock_guard<mutex> lkc{*c.origin->wait().mutex(), std::adopt_lock};

  for(size_t i=0; i<_n; ++i) { this->operator()(i) = (*c.origin)(i, c.index); }

  tock();
}

Vector Vector::Zero(size_t n, bool ready)
{
  Vector x(n);
  for(size_t i=0; i<n; ++i) { x._data[i] = 0; }
  if(ready) { x.tock(); }
  return x;
}

size_t Vector::n() const { return _n; }

double & Vector::operator()(size_t i) const
{ 
  return _data[i];
}

bool Vector::operator== (const Vector &x)
{
  WAIT_GUARD_THIS();
  WAIT_GUARD(x);

  if(x._n != _n){ return false; }

  bool result{true};
  for(size_t i=0; i<_n; ++i)
  {
    if(x._data[i] != _data[i]) { result = false; break; }  
  }

  return result;
}

Vector & Vector::operator+= (const Vector &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_plus_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Vector & Vector::operator-= (const Vector &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_sub_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Vector & Vector::operator*= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_mul_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Vector & Vector::operator/= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_div_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Scalar sven::norm(const Vector x)
{
  return sqrt(x * x);  
}

Scalar & Scalar::operator+= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_plus_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Scalar & Scalar::operator-= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_sub_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Scalar & Scalar::operator*= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_mul_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}

Scalar & Scalar::operator/= (const Scalar &x)
{
  tick();
  internal::Thunk t = [this,x]{ op_div_eq_impl(*this, x); };
  internal::RT::Q().push(t);
  return *this;
}
    
std::ostream & sven::operator<< (std::ostream &o, Vector &x)
{
  WAIT_GUARD(x);

  size_t n = x.n();
  for(size_t i=0; i<n; ++i)
  {
    o << x(i) << " ";
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
  BINOP_GUARD(a, b, ab);

  for(size_t i=0; i<a.n(); ++i) { ab(i) = a(i) + b(i); }

  ab.tock();
}

void sven::op_plus_eq_impl(Vector a, const Vector b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);

  for(size_t i=0; i<a.n(); ++i){ a(i) += b(i); }

  a.tock();
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
  BINOP_GUARD(a, b, ab);

  for(size_t i=0; i<a.n(); ++i) { ab(i) = a(i) - b(i); }

  ab.tock();
}

void sven::op_sub_eq_impl(Vector a, const Vector b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);
  
  for(size_t i=0; i<a.n(); ++i){ a(i) -= b(i); }

  a.tock();
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
  for(size_t i=0; i<n; ++i) { d += a[i*a_stride] * b[i*b_stride]; }
  return d;
}

double sven::dot(size_t n, double *a, double *b)
{
  double d{0};
  for(size_t i=0; i<n; ++i) { d += a[i] * b[i]; }
  return d;
}

void sven::op_dot_impl(const Vector a, const Vector b, Scalar ab)
{
  BINOP_GUARD(a, b, ab);

  ab() = dot(a.n(), &a(0), &b(0));

  ab.tock();
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
  BINOP_GUARD(a, b, ab);

  for(size_t i=0; i<a.n(); ++i) { ab(i) = a(i) / b(); }

  ab.tock();
}

void sven::op_div_eq_impl(Vector a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);
  
  for(size_t i=0; i<a.n(); ++i) { a(i) /= b(); }

  a.tock();
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
  BINOP_GUARD(a, b, ab);

  for(size_t i=0; i<a.n(); ++i) { ab(i) = a(i) * b(); }

  ab.tock();
}

void sven::op_mul_eq_impl(Vector a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);
  
  for(size_t i=0; i<a.n(); ++i) { a(i) *= b(); }

  a.tock();
}

//~= Scalar ~=-----------------------------------------------------------------
Scalar::Scalar() : Object(alloc<double>(1))
{
  tick();
} 

Scalar::Scalar(double value) : Scalar() 
{ 
  *_data = value; 
  tock();
}

Scalar sven::sqrt(const Scalar x)
{
  WAIT_GUARD(x);
  Scalar s(::sqrt(x()));
  return s;
}

double & Scalar::operator()() const
{ 
  return *_data;
}

bool Scalar::operator==(const Scalar &s)
{
  WAIT_GUARD_THIS();
  WAIT_GUARD(s);

  return this->operator()() == s();
}

Scalar sven::operator+ (const Scalar &a, const Scalar &b)
{
  Scalar ab;
  internal::Thunk t = [a,b,ab](){ op_plus_impl(a,b,ab); };
  internal::RT::Q().push(t);
  return ab;
}

void sven::op_plus_impl(const Scalar a, const Scalar b, Scalar ab)
{
  BINOP_GUARD(a, b, ab);

  ab() = a() + b();

  ab.tock();
}
void sven::op_plus_eq_impl(Scalar a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);

  a() += b();

  a.tock();
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

  BINOP_GUARD(a, b, ab);

  ab() = a() - b();

  ab.tock();
}
void sven::op_sub_eq_impl(Scalar a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);

  a() -= b();

  a.tock();
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
  BINOP_GUARD(a, b, ab);

  ab() = a() * b();

  ab.tock();
}
void sven::op_mul_eq_impl(Scalar a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);

  a() *= b();

  a.tock();
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
  BINOP_GUARD(a, b, ab);

  ab() = a() / b();

  ab.tock();
}
void sven::op_div_eq_impl(Scalar a, const Scalar b)
{
  MOD_GUARD(a);
  WAIT_GUARD(b);

  a() /= b();

  a.tock();
}

//~= Matrix ~=-----------------------------------------------------------------

Matrix::Matrix(size_t m, size_t n)
  : Object(alloc<double>(m*n)),
    _m{m}, _n{n}
{
  tick();
}

Matrix::Matrix(size_t m, size_t n, vector<double> values)
  : Matrix(m, n)
{
  if(values.size() != m*n)
  {
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  for(size_t i=0; i<m*n; ++i) { _data[i] = values[i]; }
  tock();
}

Matrix Matrix::Zero(size_t m, size_t n)
{
  Matrix A{m,n};
  for(size_t i=0; i<m*n; ++i)
  {
    A._data[i] = 0;
  }
  A.tock();
  return A;
}

Matrix Matrix::Identity(size_t m, size_t n)
{
  Matrix A = Matrix::Zero(m,n);

  for(size_t i=0; i<max(m,n); ++i){ A(i, i) = 1; }

  return A;
}

bool Matrix::operator== (const Matrix &A)
{
  WAIT_GUARD_THIS();
  WAIT_GUARD(A);

  if(A.m() != m() || A.n() != n()){ return false; }

  bool result{true};
  for(size_t i=0; i<_n*_m; ++i)
  {
    if(A._data[i] != _data[i]){ result = false; break; }
  }

  return result;
}

std::ostream & sven::operator<< (std::ostream &o, Matrix &A)
{
  WAIT_GUARD(A);

  o << setprecision(3) << fixed;
  for(size_t i=0; i<A.m(); ++i)
  {
    for(size_t j=0; j<A.n(); ++j) { o << A(i, j) << " "; }
    o << "\n";
  }
  return o;
}

//TODO: Eternal badness
Matrix & Matrix::operator*= (const Matrix &A)
{
  *this = *this * A;
  return *this;
}

double & Matrix::operator()(size_t i, size_t j) const
{
  return _data[i*_n + j];
}

size_t Matrix::m() const { return _m; }
size_t Matrix::n() const { return _n; }

Vector sven::operator* (const Matrix &A, const Vector &x)
{
  if(A.n() != x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector Ax(A.m());
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
  BINOP_GUARD(A, x, Ax);
  CountdownLatch cl{static_cast<int>(A.m())};

  for(size_t i=0; i<A.n(); ++i)
  {
    internal::Thunk t = [A, x, Ax, &cl, i]()
    { 
      multi_dot(A.n(), &A(i,0), &x(0), &Ax(i), cl); 
    };

    internal::RT::Q().push(t);
  }
  cl.wait();

  Ax.tock();
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
  BINOP_GUARD(A, B, AB);
  CountdownLatch cl{static_cast<int>(A.m()*B.n())};

  for(size_t i=0; i<A.m(); ++i)
  {
    for(size_t j=0; j<B.n(); ++j)
    {
      internal::Thunk t = [A, B, AB, &cl, i, j]()
      {
        multi_dot(A.n(), 
             &A(i,0), 1, 
             &B(0,j), B.n(), 
            &AB(i,j), 
            cl);
      };
      internal::RT::Q().push(t);
    }
  }
  cl.wait();

  AB.tock();
}

Column Matrix::C(size_t index)
{
  return Column(this, index);
}

ColumnRange Matrix::C(size_t begin, size_t end)
{
  return ColumnRange(this, begin, end);
}

//~=~ Column ~=~---------------------------------------------------------------

Column::Column(Matrix *origin, size_t index)
  : origin{origin}, index{index}
{}
    
Column & Column::operator= (const Vector &x)
{
  if(origin->m() < x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }
 
  origin->tick();
  
  internal::Thunk t = [this,x](){op_eq_impl(*this, x);};
  internal::RT::Q().push(t);

  return *this;
}

void sven::op_eq_impl(Column c, Vector x)
{
  WAIT_GUARD(x);
  lock_guard<mutex> lkc{*c.origin->mod_wait().mutex(), std::adopt_lock};

  size_t n = std::min(c.origin->m(), x.n());
  for(size_t i=0; i<n; ++i) { (*c.origin)(i, c.index) = x(i); }

  c.origin->tock();
}

using std::unique_ptr;
Column & Column::operator-= (const Vector &x)
{
  if(origin->m() < x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  origin->tick();

  internal::Thunk t = [this,x](){op_minus_eq_impl(*this, x);};
  internal::RT::Q().push(t);
  
  return *this;
}

void sven::op_minus_eq_impl(Column c, Vector x)
{
  WAIT_GUARD(x);
  lock_guard<mutex> lkc{*c.origin->mod_wait().mutex(), std::adopt_lock};
  
  size_t n = std::min(c.origin->m(), x.n());
  for(size_t i=0; i<n; ++i) { (*c.origin)(i, c.index) -= x(i); }

  c.origin->tock();
}

Column & Column::operator+= (const Vector &x)
{
  if(origin->m() < x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  origin->tick();

  internal::Thunk t = [this,x](){op_plus_eq_impl(*this, x);};
  internal::RT::Q().push(t);
  
  return *this;
}

void sven::op_plus_eq_impl(Column c, Vector x)
{
  WAIT_GUARD(x);
  lock_guard<mutex> lkc{*c.origin->mod_wait().mutex(), std::adopt_lock};

  for(size_t i=0; i<c.origin->m(); ++i) { (*c.origin)(i, c.index) += x(i); }

  c.origin->tock();
}

bool Column::operator== (const Vector &x)
{
  lock_guard<mutex> lkc{*origin->wait().mutex(), std::adopt_lock};
  WAIT_GUARD(x);
  
  if(origin->m() != x.n()) { return false; }

  bool result{true};
  for(size_t i=0; i<origin->m(); ++i)
  {
    if((*origin)(i, index) != x(i)){ result = false; break; }
  }

  return result;
}

//~=~ Column Range ~=~---------------------------------------------------------

ColumnRange::ColumnRange(Matrix *origin, size_t begin, size_t end)
  : origin{origin}, begin{begin}, end{end}
{}

size_t ColumnRange::m() const { return origin->m(); }
size_t ColumnRange::n() const { return end - begin + 1; }

Vector sven::operator* (const ColumnRange &C, const Vector &x)
{
  if(C.transposed)
  {
    if(C.m() != x.n())
    {
      throw runtime_error("non-conformal operation:" + CRIME_SCENE);
    }
    Vector Cx = Vector::Zero(C.n(), false);
    internal::Thunk t = [C,x,Cx](){ op_mul_impl_T(C,x,Cx); };
    internal::RT::Q().push(t);
    return Cx;
  }
  else
  {
    if(C.n() != x.n())
    {
      //throw runtime_error("non-conformal operation:" + CRIME_SCENE);
    }
    Vector Cx = Vector::Zero(C.m(), false);
    internal::Thunk t = [C,x,Cx](){ op_mul_impl(C,x,Cx); };
    internal::RT::Q().push(t);
    return Cx;
  }

}

Vector sven::operator* (const ColumnRange &C, const Column &x)
{
  lock_guard<mutex> lkc{*C.origin->wait().mutex(), adopt_lock};

  Vector cx = x;
  return C * cx;
}

void sven::op_mul_impl(const ColumnRange C, const Vector x, Vector Cx)
{
  lock_guard<mutex> lkc{*C.origin->wait().mutex(), adopt_lock};
  WAIT_GUARD(x);
  MOD_GUARD(Cx);
  CountdownLatch cl{static_cast<int>(C.m())};

  for(size_t i=0; i<C.m(); ++i)
  {
    internal::Thunk t = [C,x,Cx,&cl,i]()
    {
      multi_dot(C.n(),
          &(*C.origin)(i, C.begin),
          &x(0),
          &Cx(i),
          cl);
    };
    internal::RT::Q().push(t);
  }
  cl.wait();

  Cx.tock();
}

void sven::op_mul_impl_T(const ColumnRange C, const Vector x, Vector Cx)
{
  lock_guard<mutex> lkc{*C.origin->wait().mutex(), adopt_lock};
  WAIT_GUARD(x);
  MOD_GUARD(Cx);
  CountdownLatch cl{static_cast<int>(C.n())};

  for(size_t i=0; i<C.n(); ++i)
  {
    internal::Thunk t = [C,x,Cx,&cl,i]()
    {
      multi_dot(C.m(), 
          &(*C.origin)(0,i+C.begin), C.origin->n(),
          &x(0), 1,
          &Cx(i), cl);    
    };
    internal::RT::Q().push(t);
  }
  cl.wait();

  Cx.tock();
}

//~=~ SparseMatrix ~=~---------------------------------------------------------
  
SparseMatrixData::SparseMatrixData(size_t *r, size_t *c, double *v)
  : r{r}, c{c}, v{v}
{}

SparseMatrix::SparseMatrix(size_t m, size_t n, size_t z)
  : 
    Object(
      SparseMatrixData(
        alloc<size_t>(m),
        alloc<size_t>(m*z),
        alloc<double>(m*z))
        ),
    _m{m}, _n{n}, _z{z}
{
  tick();
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
    _data.r[i] = rs[i];
    for(size_t k=0; k<rs[i]; ++k, ++j)
    {
      if(cs[j] > n)
      {
        throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
      }
      _data.c[i*z + k] = cs[j];
      _data.v[i*z + k] = vs[j];
    }
  }

  if(j != vs.size())
  {
    throw runtime_error("unreal sparse matrix cooridnates" + CRIME_SCENE);
  }

  tock();
}

SparseMatrix SparseMatrix::Identity(size_t m, size_t n, size_t z)
{
  SparseMatrix A(m, n, z);
  for(size_t i=0; i<max(m,n); ++i)
  {
    A._data.r[i] = 1;
    A._data.c[i*z] = i;
    A._data.v[i*z] = 1;
  }
  A.tock();
  return A;
}

double & SparseMatrix::_at(size_t i, size_t j)
{
  for(size_t k=0; k<_data.r[i]; ++k)
  {
    if(_data.c[_z*i + k] == j){return _data.v[_z*i + k];}
  }

  throw runtime_error("unreal sparse matrix coordinates" + CRIME_SCENE);
}

double & SparseMatrix::operator() (size_t i, size_t j)
{
  if(i >= _m || j >= _n)
  { 
    throw runtime_error("unreal sparse matrix coordinates" + CRIME_SCENE);
  }
  return _at(i, j);
}

size_t SparseMatrix::m() const { return _m; }
size_t SparseMatrix::n() const { return _n; }
size_t SparseMatrix::z() const { return _z; }

size_t* SparseMatrix::r() const { return _data.r; }
size_t* SparseMatrix::c() const { return _data.c; }
double* SparseMatrix::v() const { return _data.v; }

Vector sven::operator* (const SparseMatrix &A, const Vector &x)
{
  if(A.n() != x.n())
  { 
    throw runtime_error("non-conformal operation:" + CRIME_SCENE); 
  }

  Vector Ax(A.m());
  internal::Thunk t = [A,x,Ax](){ op_mul_impl(A,x,Ax); };
  internal::RT::Q().push(t);
  return Ax;
}

Vector sven::operator* (const SparseMatrix &A, const Column &x)
{
  Vector Ax(A.m());
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
  for(size_t i=0; i<rz; ++i) { d += A_v[i] * x[A_c[i]]; }
  *Ax = d;
  --cl;
}

void sven::op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax)
{
  BINOP_GUARD(A, x, Ax);
  CountdownLatch cl{static_cast<int>(A.m())};

  for(size_t i=0; i<A.m(); ++i)
  {
    internal::Thunk t = [A,x,Ax,&cl,i]()
    {
      multi_sparse_dot(A.r()[i], 
          &A.c()[i*A.z()], &A.v()[i*A.z()], &x(0), &Ax(i), 
          cl);
    };
    internal::RT::Q().push(t);
  }
  cl.wait();

  Ax.tock();
}

void sven::op_mul_impl(const SparseMatrix A, const Column cx, Vector Ax)
{
  WAIT_GUARD(A);
  lock_guard<mutex> lkcx{*cx.origin->wait().mutex(), adopt_lock};
  MOD_GUARD(Ax);
  CountdownLatch cl{static_cast<int>(A.n())};

  Vector x = cx;
  for(size_t i=0; i<A.m(); ++i)
  {
    internal::Thunk t = [A,x,Ax,&cl,i]()
    {
      multi_sparse_dot(A.r()[i], 
          &A.c()[i*A.z()], &A.v()[i*A.z()], &x(0), &Ax(i), 
          cl);
    };
    internal::RT::Q().push(t);
  }
  cl.wait();

  Ax.tock();
}
