#ifndef SVEN_RUNTIME_HXX
#define SVEN_RUNTIME_HXX

#include <functional>
#include <tbb/concurrent_queue.h>
#include <thread>
#include <vector>

namespace sven { namespace internal
{

using Thunk = std::function<void()>;
using Queue = tbb::concurrent_bounded_queue<Thunk>;

class RT
{
  public:
    RT();

    static RT& get()
    {
      static RT instance;
      return instance;
    }

    static Queue & Q(){ return get()._Q; }
    std::vector<std::thread*> threads;
    void thread_func();

    Queue _Q{};
};

}}

#endif
