/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 15 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/Runtime.hxx"

using namespace sven;
using namespace sven::internal;
using std::thread;

RT::RT()
{
  size_t n = std::thread::hardware_concurrency()*4;
  threads.reserve(n);
  for(size_t i=0; i<n; ++i)
  {
    thread *t = new thread(&RT::thread_func, this);
    threads.push_back(t);
    t->detach();
  }
}

void RT::thread_func()
{
  while(47)
  {
    Thunk t;
    _Q.pop(t);
    ++active_jobs;
    t();
    --active_jobs;
  }
}
