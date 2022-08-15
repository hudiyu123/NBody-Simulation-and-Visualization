#ifndef CPP_PROJECT_INCLUDE_NBODY_CONCURRENCY_QUEUE_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_CONCURRENCY_QUEUE_HPP_

#include <concepts>
#include <condition_variable>
#include <queue>
#include <mutex>

namespace NBody::Concurrency {

/// \brief Concurrent bounded FIFO queue class
template<std::semiregular T>
class Queue {
 public:
  /// \brief The type of each of the elements stored in the queue.
  using ValueType = T;

  /// \brief A type for the status of a queue operation.
  enum class Status {
    /// \brief operation successful
    Success = 0,
    /// \brief queue is closed
    Closed,
  };

 private:
  std::queue<T> queue_;
  bool is_closed_;
  std::mutex m_;
  std::condition_variable empty_;

 public:

  /// \brief Constructs a queue.
  /// \details The queue is marked as open (i.e., not closed).
  Queue() : queue_{}, is_closed_{}, m_{}, empty_{} {}

  /// \brief A queue is not movable or copyable.
  Queue(const Queue&) = delete;
  /// \brief A queue is not movable or copyable.
  Queue& operator=(const Queue&) = delete;
  /// \brief A queue is not movable or copyable.
  Queue(Queue&&) = delete;
  /// \brief A queue is not movable or copyable.
  Queue& operator=(Queue&&) = delete;

  /// \brief Destroys the queue after closing the queue (if not already closed)
  /// and clearing the queue (if not already empty).
  ~Queue() {
    close();
    clear();
  }

  /// \brief Inserts the value x at the end of the queue, blocking if necessary.
  /// \details If the value x is successfully inserted on the queue, the
  /// function returns status::success.\n
  /// \n
  /// If the value x cannot be inserted on the queue (due to the queue being
  /// closed), the function returns with a return value of status::closed.
  /// \note This function is thread safe.
  /// \note The rvalue reference parameter is intentional and implies that the
  /// push function is permitted to change the value of x (e.g., by moving from
  /// x).
  /// \param x The value to be pushed
  /// \return The status of the push operation
  Status push(ValueType&& x) {
    auto lock = std::unique_lock{m_};
    if (is_closed()) {
      return Status::Closed;
    } else {
      queue_.push(x);
      empty_.notify_one();
      return Status::Success;
    }
  }

  /// \brief Removes the value from the front of the queue and places it in x,
  /// blocking if necessary.
  /// \details If the queue is empty and not closed, the thread is blocked
  /// until: 1) a value can be removed from the queue; or 2) the queue is
  /// closed.\n
  /// \n
  /// If the queue is closed, the function does not block and either returns
  /// status::closed or status::success, depending on whether a value can be
  /// successfully removed from the queue.\n
  /// \n
  /// If a value is successfully removed from the queue, the value is placed in
  /// x and the function returns status::success. If a value cannot be
  /// successfully removed from the queue (due to the queue being both empty and
  /// closed), the function returns status::closed.
  /// \note This function is thread safe.
  /// \param x The reference that stores the popped value
  /// \return The status of the pop operation
  Status pop(ValueType& x) {
    auto lock = std::unique_lock{m_};
    if (is_closed()) {
      if (is_empty()) {
        return Status::Closed;
      } else {
        x = queue_.front();
        queue_.pop();
        return Status::Success;
      }
    } else {
      if (is_empty()) {
        empty_.wait(lock, [this]() {
          return !is_empty() || is_closed();
        });

        if (is_closed()) {
          if (is_empty()) {
            return Status::Closed;
          } else {
            x = queue_.front();
            queue_.pop();
            return Status::Success;
          }
        } else {
          x = queue_.front();
          queue_.pop();
          return Status::Success;
        }
      } else {
        x = queue_.front();
        queue_.pop();
        return Status::Success;
      }
    }
  }

  /// \brief Closes the queue.
  /// \details The queue is placed in the closed state.\n
  /// The closed state prevents more items from being inserted on the queue, but
  /// it does not clear the items that are already on the queue.\n
  /// Invoking this function on a closed queue has no effect.\n
  /// \note This function is thread safe.
  void close() {
    auto lock = std::unique_lock{m_};
    if (!is_closed()) {
      is_closed_ = true;
      empty_.notify_all();
    }
  }

  /// \brief Clears the queue.
  /// \details All elements on the queue are discarded.\n
  /// \note This function is thread safe.
  void clear() {
    auto lock = std::unique_lock{m_};
    while (!queue_.empty()) {
      queue_.pop();
    }
  }

  /// \brief Returns if the queue is currently empty.
  /// \note This function is not thread safe.
  /// \return true if the queue is empty, false otherwise
  [[nodiscard]] bool is_empty() const {
    return queue_.empty();
  }

  /// \brief Returns if the queue is closed (i.e., in the closed state).
  /// \note This function is not thread safe.
  /// \return true if the queue is closed, false otherwise
  [[nodiscard]] bool is_closed() const {
    return is_closed_;
  }
};

} // namespace NBody::Concurrency

#endif //CPP_PROJECT_INCLUDE_NBODY_CONCURRENCY_QUEUE_HPP_
