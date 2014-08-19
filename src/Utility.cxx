/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 15 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/Utility.hxx"

using namespace sven;
using std::unique_lock;
using std::mutex;


CountdownLatch::CountdownLatch(int size)
  :_cnt{size}
{ }

void CountdownLatch::wait()
{
  unique_lock<mutex> lk{*_mtx};
  if(_cnt > 0) {
    _cnd->wait(lk);
  }
  lk.unlock();
}

void CountdownLatch::set(int count){ _cnt = count; }

void CountdownLatch::operator--()
{
  --_cnt;
  if(_cnt <= 0){_cnd->notify_all();}
  if(_cnt < 0){ _cnt = 0; }
}

void CountdownLatch::operator--(int)
{
  --_cnt;
  if(_cnt <= 0){_cnd->notify_all();}
  if(_cnt < 0){ _cnt = 0; }
}

void CountdownLatch::operator++()
{
  ++_cnt;
}

void CountdownLatch::operator++(int)
{
  ++_cnt;
}

int CountdownLatch::operator()() const
{
  return _cnt;
}

#ifdef DEBUG
void sven::log_recycle(const char *fn, size_t linno)
{
  std::stringstream ss;
  ss
    << std::hex << std::this_thread::get_id() << std::dec << "  "
    << BLU << fn << ":" GRY << linno << NRM
    << GRN << " ReCyClEd" << NRM;
  LOGG(ss.str());
}
#endif
