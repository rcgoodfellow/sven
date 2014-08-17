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

class Vector
{
  public:

    explicit Vector(size_t n);
    explicit Vector(std::initializer_list<double>);
    Vector(const Column &c);
    static Vector Zero(size_t n);

    size_t n() const;
    double & operator()(size_t i);
    ObjectState state() const;
    bool operator==(const Vector &x);
   
    Vector & operator+=(const Vector &x);
    Vector & operator-=(const Vector &x);
    Vector & operator*=(const Scalar &x);
    Vector & operator/=(const Scalar &x);


  private:
    size_t _n;
    double *_;
    std::shared_ptr<ObjectState> _state;
    std::shared_ptr<std::mutex> _mtx;
    std::shared_ptr<std::condition_variable> _cnd;


    friend void op_plus_impl(const Vector a, const Vector b, Vector ab);
    friend void op_sub_impl(const Vector a, const Vector b, Vector ab);
    friend void op_dot_impl(const Vector a, const Vector b, Scalar ab);
    friend void op_div_impl(const Vector a, const Scalar b, Vector ab);
    friend void op_mul_impl(const Vector a, const Scalar b, Vector ab);
    
    friend void op_mul_impl(const Matrix A, const Vector x, Vector mx);
    
    friend void op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax);
    
    friend void op_mul_impl(const ColumnRange C, const Vector x, Vector Cx);
    friend void op_mul_impl_T(const ColumnRange C, const Vector x, Vector Cx);

    friend Scalar norm(const Vector x);

    friend class Column;
    friend class OperandStasis<Vector>;
    friend class OperandStasis<Vector,Matrix>;
    friend class OperandStasis<Matrix,Vector>;
    friend class OperandStasis<Scalar,Vector>;
    friend class OperandStasis<Vector,Scalar>;
    friend class OperandStasis<Vector,SparseMatrix>;
    friend class OperandStasis<SparseMatrix,Vector>;

    friend std::ostream & operator<< (std::ostream &, Vector &x);
};
    
Scalar norm(const Vector x);
    
std::ostream & operator<< (std::ostream &, Vector &x);

Vector operator+ (const Vector &a, const Vector &b);
void op_plus_impl(const Vector a, const Vector b, Vector ab);

Vector operator- (const Vector &a, const Vector &b);
void op_sub_impl(const Vector a, const Vector b, Vector ab);

double dot(size_t n, double *a, size_t a_stride, double *b, size_t b_stride);
double dot(size_t n, double *a, double *b);
Scalar operator* (const Vector &a, const Vector &b);
void op_dot_impl(const Vector a, const Vector b, Scalar ab);

Vector operator/ (const Vector &a, const Scalar &b);
void op_div_impl(const Vector a, const Scalar b, Vector ab);

Vector operator* (const Vector &a, const Scalar &b);
void op_mul_impl(const Vector a, const Scalar b, Vector ab);

void multi_dot(size_t n, 
    double *a, size_t a_stride,
    double *b, size_t b_stride,
    double *ab,
    CountdownLatch &cl);

void multi_dot(size_t n, double *a, double *b, double *ab,
    CountdownLatch &cl);

class Scalar
{
  public:
    Scalar();
    explicit Scalar(double value);

    double & operator()();
    bool operator==(const Scalar &s);
    ObjectState state() const;

    Scalar & operator+= (const Scalar&);
    Scalar & operator-= (const Scalar&);
    Scalar & operator*= (const Scalar&);
    Scalar & operator/= (const Scalar&);

  private:
    double *_;
    ObjectState *_state{new ObjectState};
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};

    friend void op_div_impl(const Vector a, const Scalar b, Vector ab);
    friend void op_mul_impl(const Vector a, const Scalar b, Vector ab);
    friend void op_dot_impl(const Vector a, const Vector b, Scalar ab);

    friend void op_plus_impl(const Scalar a, const Scalar b, Scalar ab);
    friend void op_sub_impl (const Scalar a, const Scalar b, Scalar ab);
    friend void op_mul_impl(const Scalar a, const Scalar b, Scalar ab);
    friend void op_div_impl(const Scalar a, const Scalar b, Scalar ab);
    
    friend class OperandStasis<Scalar>;
    friend class OperandStasis<Scalar,Vector>;
    friend class OperandStasis<Vector,Scalar>;

    friend Scalar sqrt(const Scalar);
    
};

Scalar sqrt(const Scalar);

Scalar operator+ (const Scalar &a, const Scalar &b);
void op_plus_impl(const Scalar a, const Scalar b, Scalar ab);

Scalar operator- (const Scalar &a, const Scalar &b);
void op_sub_impl (const Scalar a, const Scalar b, Scalar ab);

Scalar operator* (const Scalar &a, const Scalar &b);
void op_mul_impl(const Scalar a, const Scalar b, Scalar ab);

Scalar operator/ (const Scalar &a, const Scalar &b);
void op_div_impl (const Scalar a, const Scalar b, Scalar ab);

class Column
{
  public:
    Column() = delete;
    Column & operator= (const Vector &);
    bool operator== (const Vector &);

  private:
    Column(Matrix *origin, size_t index);
    Matrix *_origin;
    size_t _index;

    friend class Matrix;
    friend class Vector;
};

class Matrix
{
  public:
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, std::vector<double> values);
    static Matrix Zero(size_t m, size_t n);
    static Matrix Identity(size_t m, size_t n);
    bool operator==(const Matrix &A);
    
    Matrix & operator*= (const Matrix &);

    double & operator()(size_t i, size_t j);
    Column C(size_t index);
    ColumnRange C(size_t begin, size_t end);

    size_t m() const, n() const;
    ObjectState state() const;

  private:
    size_t _m, _n;
    double *_;
    ObjectState *_state{new ObjectState};
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};

    friend void op_mul_impl(const Matrix A, const Vector x, Vector Ax);
    friend void op_mul_impl(const Matrix A, const Matrix B, Matrix AB);

    friend void op_mul_impl(const ColumnRange C, const Vector x, Vector Cx);
    friend void op_mul_impl_T(const ColumnRange C, const Vector x, Vector Cx);

    friend class Vector;
    friend class Column;
    friend class OperandStasis<Matrix>;
    friend class OperandStasis<Matrix,Vector>;
    friend class OperandStasis<Vector,Matrix>;

    friend std::ostream & operator<< (std::ostream &, Matrix &x);
};

std::ostream & operator<< (std::ostream &, Matrix &x);

Vector operator* (const Matrix &A, const Vector &x);
void op_mul_impl(const Matrix A, const Vector x, Vector Ax);

Matrix operator* (const Matrix &A, const Matrix &B);
void op_mul_impl(const Matrix A, const Matrix B, Matrix AB);

class ColumnRange
{
  public:
    ColumnRange() = delete;
    ColumnRange(Matrix *origin, size_t begin, size_t end);

    size_t m() const, n() const;

    bool transposed{false};

  private:
    Matrix *_origin;
    size_t _begin, _end;

    friend void op_mul_impl(const ColumnRange C, const Vector x, Vector Cx);
    friend void op_mul_impl_T(const ColumnRange C, const Vector x, Vector Cx);
};

Vector operator* (const ColumnRange &C, const Vector &x);
Vector operator* (const ColumnRange &C, const Column &x);
void op_mul_impl(const ColumnRange C, const Vector x, Vector Cx);
void op_mul_impl_T(const ColumnRange C, const Vector x, Vector Cx);

class SparseMatrix
{
  public:
    SparseMatrix(size_t m, size_t n, size_t z);
    SparseMatrix(size_t m, size_t n, size_t z, 
        std::vector<size_t> rs, std::vector<size_t> cs, 
        std::vector<double> vs);
    static SparseMatrix Identity(size_t m, size_t n, size_t z);

    double & operator()(size_t i, size_t j);

    size_t m() const, n() const, z() const;
    ObjectState state() const;

  private:
    size_t _m, _n, _z;
    size_t *_r, *_c;
    double *_v;
    ObjectState *_state{new ObjectState};
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};

    double & _at(size_t i, size_t j);

    friend void op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax);
    
    friend class OperandStasis<SparseMatrix>;
    friend class OperandStasis<Vector,SparseMatrix>;
    friend class OperandStasis<SparseMatrix,Vector>;
    
};

void multi_sparse_dot(size_t z,
    size_t *A_c, double *A_v,
    double *x,
    double *Ax,
    CountdownLatch &cl);

Vector operator* (const SparseMatrix &A, const Vector &x);
void op_mul_impl(const SparseMatrix A, const Vector x, Vector Ax);

Vector operator* (const SparseMatrix &A, const Column &x);

}
#endif
