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

class Scalar;

class Vector
{
  public:

    explicit Vector(size_t n);
    explicit Vector(std::initializer_list<double>);
    static Vector Zero(size_t n);

    size_t n() const;
    double & operator()(size_t i);
    ObjectState state() const;

  private:
    size_t _n;
    double *_;
    ObjectState *_state{new ObjectState};
    std::mutex *_mtx{new std::mutex};
    std::condition_variable *_cnd{new std::condition_variable};

    friend void op_plus_impl(const Vector a, const Vector b, Vector ab);
    friend void op_dot_impl(const Vector a, const Vector b, Scalar ab);
};


Vector operator+ (const Vector &a, const Vector &b);
void op_plus_impl(const Vector a, const Vector b, Vector ab);

Scalar operator* (const Vector &a, const Vector &b);
void op_dot_impl(const Vector a, const Vector b, Scalar ab);

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
    std::mutex *_mtx{new std::mutex};
    std::condition_variable *_cnd{new std::condition_variable};

    friend void op_dot_impl(const Vector a, const Vector b, Scalar ab);
};

class Matrix
{
  public:
    Matrix(size_t m, size_t n);
    static Matrix Zero(size_t m, size_t n);
    static Matrix Identity(size_t m, size_t n);

    size_t m(), n();

  private:
    size_t _m, _n;
    double *_;
};

}
#endif
