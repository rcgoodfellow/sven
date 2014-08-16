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

namespace sven
{

enum class ObjectState{Materializing, SolidState, Vapor};

class Vector;
class Scalar;
class Matrix;

class Vector
{
  public:

    explicit Vector(size_t n);
    explicit Vector(std::initializer_list<double>);
    static Vector Zero(size_t n);

    //Vector(const Vector &);
    //Vector & operator=(const Vector &);

    size_t n() const;
    double & operator()(size_t i);
    ObjectState state() const;
    bool operator==(const Vector &x);

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
};

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
    ObjectState state() const;

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
};

Scalar operator+ (const Scalar &a, const Scalar &b);
void op_plus_impl(const Scalar a, const Scalar b, Scalar ab);

Scalar operator- (const Scalar &a, const Scalar &b);
void op_sub_impl (const Scalar a, const Scalar b, Scalar ab);

Scalar operator* (const Scalar &a, const Scalar &b);
void op_mul_impl(const Scalar a, const Scalar b, Scalar ab);

Scalar operator/ (const Scalar &a, const Scalar &b);
void op_div_impl (const Scalar a, const Scalar b, Scalar ab);

class Matrix
{
  public:
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, std::vector<double> values);
    static Matrix Zero(size_t m, size_t n);
    static Matrix Identity(size_t m, size_t n);
    bool operator==(const Matrix &A);

    size_t m() const, n() const;

  private:
    size_t _m, _n;
    double *_;
    ObjectState *_state{new ObjectState};
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};

    friend void op_mul_impl(const Matrix A, const Vector x, Vector Ax);
    friend void op_mul_impl(const Matrix A, const Matrix B, Matrix AB);
};

Vector operator* (const Matrix &A, const Vector &x);
void op_mul_impl(const Matrix A, const Vector x, Vector Ax);

Matrix operator* (const Matrix &A, const Matrix &B);
void op_mul_impl(const Matrix A, const Matrix B, Matrix AB);

}
#endif
