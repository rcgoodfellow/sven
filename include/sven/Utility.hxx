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

#define SVEN_DEFAULT_ALIGN 64

#define CRIME_SCENE string(__FILE__) + string(":") + to_string(__LINE__) + string(":") + string(__func__)

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

  private:
    std::mutex *_mtx;
    std::condition_variable *_cnd;
};

}

#endif
