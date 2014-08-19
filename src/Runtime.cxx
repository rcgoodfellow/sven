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
using std::mutex;

mutex RT::print_mtx{};

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
  //_Q.set_capacity(n/2);
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

size_t RT::next_id()
{
#ifdef DEBUG
  name_map[obj_counter] = std::to_string(obj_counter);
#endif

  return obj_counter++;
}
