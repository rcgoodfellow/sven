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
  _cnd->wait(lk);
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
  this->operator--();
}

void CountdownLatch::operator++()
{
  _cnt++;
}

void CountdownLatch::operator++(int)
{
  _cnt++;
}
