#ifndef HCI_TIMER_H_
#define HCI_TIMER_H_

#include "std.h"

#include "parallel.h"

class Time {
 private:
  static std::vector<int>& get_index() {
    static std::vector<int> index(1, 0);
    return index;
  };

  // Event index string to be printed with each start and end information.
  static std::string get_index_string() {
    const int INDEX_IGNORE = 2;
    auto& index = Time::get_index();
    if (INDEX_IGNORE >= index.size() - 1) return "";
    std::stringstream ss;
    for (std::size_t i = INDEX_IGNORE; i < index.size() - 1; i++) ss << index[i] << '.';
    return ss.str() + " ";
  }

  static std::chrono::high_resolution_clock::time_point& get_timer(const std::string& event) {
    using namespace std::chrono;
    static std::unordered_map<std::string, high_resolution_clock::time_point> times;
    return times[event];
  }

 public:
  static void init() {
    using namespace std::chrono;

    if (Parallel::get_id() != 0) return;
    Time::get_timer("init") = high_resolution_clock::now();
  }

  static void start(const std::string& event) {
    using namespace std::chrono;

    Parallel::barrier();
    if (Parallel::get_id() != 0) return;

    const auto now = high_resolution_clock::now();
    const auto init_time = Time::get_timer("init");
    Time::get_timer(event) = now;
    const duration<double> since_init = duration_cast<duration<double>>(now - init_time);
    auto& index = Time::get_index();
    index.back()++;
    index.push_back(0);
    printf(
        "\n%sBEGIN %s [%.3f/0.000]\n",
        Time::get_index_string().c_str(),
        event.c_str(),
        since_init.count());
  };

  static void end(const std::string& event) {
    using namespace std::chrono;

    Parallel::barrier();
    if (Parallel::get_id() != 0) return;

    const auto now = high_resolution_clock::now();
    const auto start_time = Time::get_timer(event);
    const auto init_time = Time::get_timer("init");
    const duration<double> since_start = duration_cast<duration<double>>(now - start_time);
    const duration<double> since_init = duration_cast<duration<double>>(now - init_time);
    printf(
        "%sEND %s [%.3f/%.3f] \n",
        Time::get_index_string().c_str(),
        event.c_str(),
        since_init.count(),
        since_start.count());
    Time::get_index().pop_back();
  }
};

#endif