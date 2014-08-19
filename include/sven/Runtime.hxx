#ifndef SVEN_RUNTIME_HXX
#define SVEN_RUNTIME_HXX

/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * This file contains the core runtime class function and template declarations
 *
 * 15 August 2014
 * ~ ry
 *
 *****************************************************************************/

#include <functional>
#include <tbb/concurrent_queue.h>
#include <thread>
#include <vector>
#include <iostream>
#include <string>
#include <unordered_map>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

namespace sven { namespace internal
{

using Thunk = std::function<void()>;
using Queue = tbb::concurrent_bounded_queue<Thunk>;

#ifdef DEBUG
#define LOGG(_MSG_) sven::internal::RT::log(_MSG_);
#else
#define LOGG(_MSG_) ;
#endif

class RT
{
  public:
    RT();

    size_t obj_counter{0};

    size_t next_id();

#ifdef DEBUG
    std::unordered_map<size_t, std::string> name_map;

    template<class T>
    static void map_name(T obj, std::string name)
    {
      get().name_map[obj.id()] = name;
    }

    template<class T>
    static std::string get_name(T obj)
    {

      try{ return get().name_map.at(obj.id()); }
      catch(...) 
      { 
        //void *arr[10];
        //size_t size;
    
        //size = backtrace(arr, 10);
        //backtrace_symbols_fd(arr, size, STDERR_FILENO);
        //exit(1);
        throw;
      }
    }
#endif

    static RT& get()
    {
      static RT instance;
      return instance;
    }

    static Queue & Q(){ return get()._Q; }
    std::vector<std::thread*> threads;
    void thread_func();

    Queue _Q{};
    std::atomic<int> active_jobs{0};

    static std::mutex print_mtx;

    static void log(std::string s)
    {
      print_mtx.lock();
      std::cout << s << std::endl;
      print_mtx.unlock();
    }
};

}}

#endif
