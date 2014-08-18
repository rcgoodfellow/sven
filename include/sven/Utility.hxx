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

#include <iostream>
#include <string>

namespace sven {

//Allocate a typed, aligned chunk of memory without havint to do the typecast
//and sizeof math bit each time
template<class T>
T* alloc(size_t sz)
{
  return (T*)_mm_malloc(sizeof(T)*sz, SVEN_DEFAULT_ALIGN);
}

#define WAIT_GUARD(__X__) \
  std::lock_guard<std::mutex> lk##__X__{*__X__.wait().mutex(), std::adopt_lock};

#define WAIT_GUARD_THIS() \
  std::lock_guard<std::mutex> lkme{*this->wait().mutex(), std::adopt_lock};

#define MOD_GUARD(__X__) \
  std::lock_guard<std::mutex> lk##__X__{*__X__.mod_wait().mutex(), \
                                                          std::adopt_lock};

#define HOLD_GUARD(__X__) \
  std::lock_guard<std::mutex> lk##__X__{*__X__.hold().mutex(), std::adopt_lock};

#define BINOP_GUARD(__A__, __B__, __AB__) \
  WAIT_GUARD(__A__); \
  WAIT_GUARD(__B__); \
  MOD_GUARD(__AB__);


template<class T>
class Object
{
  public:
    size_t tick()
    {
      return ++_scheduled_version;
    }
    
    size_t tock()
    {
      ++(*_current_version);
      _cnd->notify_all();
      return *_current_version;
    }

    bool up_to_date() const
    {
      return *_current_version == _scheduled_version;
    }

    std::unique_lock<std::mutex> wait() const
    {
      std::unique_lock<std::mutex> lk{*_mtx};
      if(!up_to_date())
      {
        _cnd->wait(lk, 
            [this](){return *_current_version == _scheduled_version;});
      }
      return lk;
    }

    std::unique_lock<std::mutex> mod_wait() const
    {
      std::unique_lock<std::mutex> lk{*_mtx};
      if(*_current_version <= (_scheduled_version - 1))
      {
        _cnd->wait(lk, 
            [this](){return *_current_version == (_scheduled_version - 1);});
      }
      return lk;
    }

    std::unique_lock<std::mutex> hold() const
    {
      return std::unique_lock<std::mutex>{*_mtx};
    }

    Object(T data) : _data{data}
    {
      *_current_version = 0;
    }

  protected:
    T _data;
    size_t *_current_version{new size_t};
    size_t _scheduled_version{0};

    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};
};

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

template<class E> 
E& T(E &e){ e.transposed = true; return e; }

template<class E> 
E T(E e){ e.transposed = true; return e; }

}

#endif
