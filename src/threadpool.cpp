/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <kmdiff/threadpool.hpp>

namespace kmdiff {

  ThreadPool::ThreadPool(size_type threads)
  {
    if (threads < _n) _n = threads;
    for (int i = 0; i < _n; i++)
    {
      _pool.push_back(std::thread(&ThreadPool::worker, this, i));
    }
  }

  ThreadPool::~ThreadPool()
  {
    {
      std::unique_lock<std::mutex> lock(_queue_mutex);
      _stop = true;
    }
    _condition.notify_all();
    for (std::thread& t : _pool)
      if (t.joinable()) t.join();
  }

  void ThreadPool::join_all()
  {
    {
      std::unique_lock<std::mutex> lock(_queue_mutex);
      _stop = true;
    }
    _condition.notify_all();
    for (std::thread& t : _pool) t.join();
  }

  void ThreadPool::join(int i)
  {
    if (_pool[i].joinable()) _pool[i].join();
  }

  void ThreadPool::worker(int i)
  {
    while (true)
    {
      std::function<void(int)> task;
      {
        std::unique_lock<std::mutex> lock(this->_queue_mutex);
        this->_condition.wait(lock, [this] { return this->_stop || !this->_queue.empty(); });
        if (this->_stop && this->_queue.empty()) return;
        task = std::move(this->_queue.front());
        this->_queue.pop();
      }
      task(i);
    }
  }

} // end of namespace kmdiff

