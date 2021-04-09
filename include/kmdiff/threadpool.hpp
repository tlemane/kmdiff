/*****************************************************************************
 *   kmtricks-sv
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

#pragma once

// std
#include <functional>
#include <future>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace kmdiff
{
class ThreadPool
{
  using size_type = std::result_of<decltype (&std::thread::hardware_concurrency)()>::type;

 public:
  ThreadPool(size_type threads);

  ThreadPool() = delete;
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;
  ThreadPool(ThreadPool&&) = delete;
  ThreadPool& operator=(ThreadPool&&) = delete;

  ~ThreadPool();

  void join_all();

  void join(int i);

  template <typename Callable>
  void add_task(Callable&& f)
  {
    auto task = std::make_shared<std::packaged_task<void(int)>>(std::forward<Callable>(f));
    {
      std::unique_lock<std::mutex> lock(_queue_mutex);
      if (_stop) throw std::runtime_error("Push on stopped Pool.");
      _queue.emplace([task](int thread_id) { (*task)(thread_id); });
    }
    _condition.notify_one();
  }

 private:
  void worker(int i);

 private:
  size_type _n{std::thread::hardware_concurrency()};
  std::vector<std::thread> _pool;
  std::queue<std::function<void(int)>> _queue;
  std::mutex _queue_mutex;
  std::condition_variable _condition;
  bool _stop{false};
};

};  // namespace kmdiff