#ifndef SVEN_COREMATH_HXX
#define SVEN_COREMATH_HXX

/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * The core math header contains the core class and function declarations for
 * doing math with Sven
 *
 * 14 August 2014
 * ~ ry
 *
 *****************************************************************************/

#include "Utility.hxx"
#include "Runtime.hxx"
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace sven
{

class Vector;
class Scalar;
class Matrix;
class SparseMatrix;
class Column;
class ColumnRange;

class Vector : public Object<double*>
{
  public:

    explicit Vector(size_t n);
    explicit Vector(std::initializer_list<double>);
    Vector(const Column &c);
    static Vector Zero(size_t n, bool ready=true);
    Vector operator! ();


    size_t n() const;
    double & operator()(size_t i) const;
    bool operator==(const Vector &x);
   
    Vector & operator+=(const Vector &x);
    Vector & operator-=(const Vector &x);
    Vector & operator*=(const Scalar &x);
    Vector & operator/=(const Scalar &x);


  private:
    size_t _n;
};
    
Scalar norm(const Vector x);
    
std::ostream & operator<< (std::ostream &, Vector &x);

void op_bang_impl(Vector a, const Vector b);

Vector operator+ (const Vector &a, const Vector &b);
void op_plus_impl(const Vector a, const Vector b, Vector ab);
void op_plus_eq_impl(Vector a, const Vector b);

Vector operator- (const Vector &a, const Vector &b);
void op_sub_impl(const Vector a, const Vector b, Vector ab);
void op_sub_eq_impl(Vector a, const Vector b);

Vector operator/ (const Vector &a, const Scalar &b);
void op_div_impl(const Vector a, const Scalar b, Vector ab);
void op_div_eq_impl(Vector a, const Scalar b);

Vector operator* (const Vector &a, const Scalar &b);
void op_mul_impl(const Vector a, const Scalar b, Vector ab);
void op_mul_eq_impl(Vector a, const Scalar b);

Scalar operator* (const Vector &a, const Vector &b);
double dot(size_t n, double *a, size_t a_stride, double *b, size_t b_stride);
double dot(size_t n, double *a, double *b);
void op_dot_impl(const Vector a, const Vector b, Scalar ab);


void multi_dot(size_t n, 
    double *a, size_t a_stride,
    double *b, size_t b_stride,
    double *ab,
    CountdownLatch &cl);

void multi_dot(size_t n, double *a, double *b, double *ab,
    CountdownLatch &cl);

class Scalar : public Object<double*>
{
  public:
    Scalar();
    explicit Scalar(double value);

    double & operator()() const;
    bool operator==(const Scalar &s);

    Scalar & operator+= (const Scalar&);
    Scalar & operator-= (const Scalar&);
    Scalar & operator*= (const Scalar&);
    Scalar & operator/= (const Scalar&);

  private:
    
};

Scalar sqrt(const Scalar);

Scalar operator+ (const Scalar &a, const Scalar &b);
void op_plus_impl(const Scalar a, const Scalar b, Scalar ab);
void op_plus_eq_impl(Scalar a, const Scalar b);

Scalar operator- (const Scalar &a, const Scalar &b);
void op_sub_impl (const Scalar a, const Scalar b, Scalar ab);
void op_sub_eq_impl(Scalar a, const Scalar b);

Scalar operator* (const Scalar &a, const Scalar &b);
void op_mul_impl(const Scalar a, const Scalar b, Scalar ab);
void op_mul_eq_impl(Scalar a, const Scalar b);

Scalar operator/ (const Scalar &a, const Scalar &b);
void op_div_impl (const Scalar a, const Scalar b, Scalar ab);
void op_div_eq_impl(Scalar a, const Scalar b);

class Column
{
  public:
    Column() = delete;
    Column & operator= (const Vector &);
    Column & operator-= (const Vector &);
    Column & operator+= (const Vector &);
    bool operator== (const Vector &);
    double & operator()(size_t i);
    
    Matrix *origin;
    size_t index;

  private:
    Column(Matrix *origin, size_t index);

    friend class Matrix;
    friend class Vector;

};

void op_eq_impl(Column c, Vector x);
void op_minus_eq_impl(Column c, Vector x);
void op_plus_eq_impl(Column c, Vector x);

class Matrix : public Object<double*>
{
  public:
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, std::vector<double> values);
    static Matrix Zero(size_t m, size_t n);
    static Matrix Identity(size_t m, size_t n);
    bool operator==(const Matrix &A);

    Matrix operator!();
    
    Matrix & operator*= (const Matrix &);

    double & operator()(size_t i, size_t j) const;
    Column C(size_t index);
    ColumnRange C(size_t begin, size_t end);

    size_t m() const, n() const;

  private:
    size_t _m, _n;
};

std::ostream & operator<< (std::ostream &, Matrix &x);

void op_bang_impl(Matrix A, const Matrix B, bool use_mod_guard = false);

Vector operator* (const Matrix &A, const Vector &x);
void op_mul_impl(const Matrix A, const Vector x, Vector Ax);

Matrix operator* (const Matrix &A, const Matrix &B);
void op_mul_impl(const Matrix A, const Matrix B, Matrix AB);
void op_mul_eq_impl(Matrix A, const Matrix B);

class ColumnRange
{
  public:
    ColumnRange() = delete;
    ColumnRange(Matrix *origin, size_t begin, size_t end);

    size_t m() const, n() const;

    bool transposed{false};
    
    Matrix *origin;
    size_t begin, end;

  private:

};

Vector operator* (const ColumnRange &C, const Vector &x);
Vector operator* (const ColumnRange &C, const Column &x);
void op_mul_impl(const ColumnRange C, size_t qcv, const Vector x, Vector Cx);
void op_mul_impl_T(const ColumnRange C, size_t qcv, const Vector x, Vector Cx);

struct SparseMatrixData
{
  size_t *r{nullptr}, *c{nullptr};
  double *v{nullptr};

  SparseMatrixData() = default;
  SparseMatrixData(size_t *r, size_t *c, double *v); 
};

class SparseMatrix : public Object<SparseMatrixData>
{
  public:
    SparseMatrix(size_t m, size_t n, size_t z);
    SparseMatrix(size_t m, size_t n, size_t z, 
        std::vector<size_t> rs, std::vector<size_t> cs, 
        std::vector<double> vs);
    static SparseMatrix Identity(size_t m, size_t n, size_t z);

    double & operator()(size_t i, size_t j);

    size_t m() const, n() const, z() const;
    
    size_t *r() const, *c() const;
    double *v() const;

  private:
    size_t _m, _n, _z;
    double & _at(size_t i, size_t j);
};

void multi_sparse_dot(size_t z,
    size_t *A_c, double *A_v,
    double *x,
    double *Ax,
    CountdownLatch &cl);

Vector operator* (const SparseMatrix &A, const Vector &x);
void op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax);

void op_mul_impl(const SparseMatrix A, const Column cx, size_t cqv, Vector Ax);

Vector operator* (const SparseMatrix &A, const Column &x);

}
#endif
