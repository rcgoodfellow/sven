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
#include <sstream>

#include "Runtime.hxx"

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

#ifdef DEBUG
template<class A>
std::string vstr(A a, size_t dc = 0, size_t ds = 0)
{
  return "[" 
    + std::to_string(a.current_version()+dc) + ":" 
    + std::to_string(a.scheduled_version()+ds) + "]";
}
#endif

#define BLU "\e[0;34m"
#define RED "\e[0;31m"
#define YEL "\e[0;33m"
#define PRP "\e[0;35m"
#define GRN "\e[0;32m"
#define CYN "\e[0;36m"
#define GRY "\e[0;30m"
#define NRM "\e[0m"

#define NAME(__X__) internal::RT::get_name(__X__)

#define MARKER(_M_) LOGG(RED _M_ NRM);

#ifdef DEBUG
template<class A, class B, class C>
void binop_log_in(const char* fn, size_t linno, const A &a, const B &b, 
    const C &c)
{
  std::stringstream ss;
  ss 
     << std::hex << std::this_thread::get_id() << std::dec << "  "
     << BLU << fn << ":" << GRY << linno << NRM << "  "
     << NAME(c) << GRY << vstr(c) << NRM
     << CYN << " <-- " << NRM
     << NAME(a) << GRY << vstr(a) << NRM <<"  "
     << NAME(b) << GRY << vstr(b) << NRM;
  LOGG(ss.str());
}

template<class A, class B>
void uop_log_in(const char* fn, size_t linno, const A &a, const B &b)
{
  std::stringstream ss;
  ss 
     << std::hex << std::this_thread::get_id() << std::dec << "  "
     << BLU << fn << ":" << GRY << linno << NRM << "  "
     << NAME(a) << GRY << vstr(a) << NRM
     << CYN << " <-- " << NRM
     << NAME(b) << GRY << vstr(b) << NRM;
  LOGG(ss.str());
}


template<class A>
void logout(const char* fn, size_t linno, const A &a)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << BLU << fn << ":" GRY << linno << GRN << "  --> " << NRM
    << NAME(a) << GRY << vstr(a) << NRM;
  LOGG(ss.str());
}

template<class A>
void wait_log(const A &a)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << CYN << "@" << NRM << NAME(a) << GRY << vstr(a) << NRM;
  LOGG(ss.str());
}

template<class A>
void unwait_log(const A &a)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << GRN << "@" << NRM << NAME(a) << GRY << vstr(a) << NRM;
  LOGG(ss.str());
}

template<class L>
void latch_log(const char * fn, size_t linno, const L &l)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << BLU << fn << ":" GRY << linno << NRM
    << CYN << " #" << GRY << l() << NRM;
  LOGG(ss.str());
}

template<class L>
void unlatch_log(const char *fn, size_t linno, const L &l)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << BLU << fn << ":" GRY << linno << NRM
    << GRN << " #" << GRY << l() << NRM;
  LOGG(ss.str());
}

void log_recycle(const char *fn, size_t linno);
#endif

#ifdef DEBUG
#define BINOP_LOG_IN(__A__, __B__, __C__) \
        binop_log_in(__func__, __LINE__, __A__, __B__, __C__);
#define UOP_LOG_IN(__A__, __B__) uop_log_in(__func__, __LINE__, __A__, __B__);
#define LOG_OUT(__A__) logout(__func__, __LINE__, __A__);
#define LATCH_LOG(__L__) latch_log(__func__, __LINE__, __L__);
#define UNLATCH_LOG(__L__) unlatch_log(__func__, __LINE__, __L__);
#define WAIT_LOG(__X__) wait_log(__X__);
#define UNWAIT_LOG(__X__) unwait_log(__X__);
#define LOG_RECYCLE() log_recycle(__func__, __LINE__);
#else
#define BINOP_LOG_IN(__A__, __B__, __C__) ;;
#define UOP_LOG_IN(__A__, __B__) ;;
#define LOG_OUT(__A__) ;;
#define LATCH_LOG(__L__) ;;
#define UNLATCH_LOG(__L__) ;;
#define WAIT_LOG(__X__) ;;
#define UNWAIT_LOG(__X__) ;;
#define LOG_RECYCLE();
#endif

#define BINOP_GUARD(__A__, __B__, __AB__) \
  WAIT_GUARD(__A__); \
  WAIT_GUARD(__B__); \
  MOD_GUARD(__AB__);

#define RECYCLE_THRESHOLD 3

template<class A, class B, class C>
bool should_recycle(A a, B b, C c)
{
  return a.vdelta() > RECYCLE_THRESHOLD || b.vdelta() > RECYCLE_THRESHOLD
    || c.vdelta() > RECYCLE_THRESHOLD;
}

template<class A, class B>
bool should_recycle(A a, B b)
{
  return a.vdelta() > RECYCLE_THRESHOLD || b.vdelta() > RECYCLE_THRESHOLD;
}

#define BINOP_RECYCLE_CHECK(__A__, __B__, __C__, __fn__)\
  if(should_recycle(__A__, __B__, __C__))\
  {\
    LOG_RECYCLE();\
    internal::Thunk t = [__A__,__B__,__C__]{ __fn__(__A__,__B__,__C__); }; \
    internal::RT::Q().push(t); \
    return; \
  }

#define BINOP_RECYCLE_CHECK_SV(__A__, __B__, __C__, __fn__, __sv__)\
  if(should_recycle(__A__, __B__, __C__))\
  {\
    LOG_RECYCLE();\
    internal::Thunk t = [__A__,__B__,__C__,__sv__]\
      { __fn__(__A__,__B__,__C__,__sv__); }; \
    internal::RT::Q().push(t); \
    return; \
  }

#define UOP_RECYCLE_CHECK(__A__, __B__,__fn__)\
  if(should_recycle(__A__, __B__))\
  {\
    LOG_RECYCLE();\
    internal::Thunk t = [__A__,__B__]{ __fn__(__A__,__B__); }; \
    internal::RT::Q().push(t); \
    return; \
  }

#define UOP_RECYCLE_CHECK_SV(__A__, __B__,__fn__,__sv__)\
  if(should_recycle(__A__, __B__))\
  {\
    LOG_RECYCLE();\
    internal::Thunk t = [__A__,__B__,__sv__]{ __fn__(__A__,__B__,__sv__); }; \
    internal::RT::Q().push(t); \
    return; \
  }


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
      return *_current_version >= _scheduled_version;
    }
    
    bool up_to_date(size_t v) const
    {
      return *_current_version >= v;
    }

    std::unique_lock<std::mutex> wait() const
    {
      WAIT_LOG(*this);
      std::unique_lock<std::mutex> lk{*_mtx};
      if(!up_to_date())
      {
        _cnd->wait(lk, 
            [this](){return *_current_version >= _scheduled_version;});
      }
      UNWAIT_LOG(*this);
      return lk;
    }
    
    std::unique_lock<std::mutex> wait(size_t v) const
    {
      WAIT_LOG(*this);
      std::unique_lock<std::mutex> lk{*_mtx};
      if(!up_to_date(v))
      {
        _cnd->wait(lk, 
            [this,v](){return *_current_version >= v;});
      }
      UNWAIT_LOG(*this);
      return lk;
    }

    std::unique_lock<std::mutex> mod_wait() const
    {
      WAIT_LOG(*this);
      std::unique_lock<std::mutex> lk{*_mtx};
      if(*_current_version <= (_scheduled_version - 1))
      {
        _cnd->wait(lk, 
            [this](){return *_current_version >= (_scheduled_version - 1);});
      }
      UNWAIT_LOG(*this);
      return lk;
    }
    std::unique_lock<std::mutex> mod_wait(size_t v) const
    {
      WAIT_LOG(*this);
      std::unique_lock<std::mutex> lk{*_mtx};
      if(*_current_version <= (v - 1))
      {
        _cnd->wait(lk, 
            [this, v](){return *_current_version >= (v - 1);});
      }
      UNWAIT_LOG(*this);
      return lk;
    }

    std::unique_lock<std::mutex> hold() const
    {
      WAIT_LOG(*this);
      return std::unique_lock<std::mutex>{*_mtx};
      UNWAIT_LOG(*this);
    }

    Object(T data) : _data{data}, _id{new size_t}
    {
      *_current_version = 0;
      *_id = internal::RT::get().next_id();
    }

    size_t scheduled_version() const { return _scheduled_version; }
    size_t current_version() const { return *_current_version; }

    long vdelta() const { return scheduled_version() - current_version(); }

    size_t id() const { return *_id; }

  protected:
    T _data;
    size_t *_current_version{new size_t};
    size_t _scheduled_version{0};

    size_t *_id;

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
    int operator()() const;

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
