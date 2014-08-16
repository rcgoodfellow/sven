#ifndef SVEN_UTILITY_HXX
#define SVEN_UTILITY_HXX

/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * This file contains usefull templates, and function declarations and classes
 * used throught Sven.
 *
 * 14 August 2014
 * ~ ry
 *
 *****************************************************************************/

#include <cstddef>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <mm_malloc.h>
#include <string>

#define SVEN_DEFAULT_ALIGN 64

#define CRIME_SCENE \
  std::string(__FILE__) + std::string(":") + \
  std::to_string(__LINE__) + \
  std::string(":") + std::string(__func__)

namespace sven {

//Allocate a typed, aligned chunk of memory without havint to do the typecast
//and sizeof math bit each time
template<class T>
T* alloc(size_t sz)
{
  return (T*)_mm_malloc(sizeof(T)*sz, SVEN_DEFAULT_ALIGN);
}

class CountdownLatch
{
  public:
    CountdownLatch(int);
    void wait();
    void set(int);
    void operator--();
    void operator--(int);
    void operator++();
    void operator++(int);

  private:
    std::atomic<int> _cnt;
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};
};

enum class ObjectState{Materializing, SolidState, Vapor};

template<class A, class B>
class OperandStasis
{
  A a;
  B b;

  std::unique_lock<std::mutex> lk_a, lk_b;

  public:
    OperandStasis(A a, B b) 
      : a{a}, b{b}, 
        lk_a{*a._mtx, std::defer_lock}, lk_b{*b._mtx, std::defer_lock}
    {
      lk_a.lock();
      if(a.state() != ObjectState::SolidState)
      {
        a._cnd->wait(lk_a);
      }

      lk_b.lock();
      if(b.state() != ObjectState::SolidState)
      {
        lk_a.unlock();
        b._cnd->wait(lk_b);
        lk_a.lock();
      }
    }

    ~OperandStasis()
    {
      lk_a.unlock();
      lk_b.unlock();
    }
  
};


}

#endif
