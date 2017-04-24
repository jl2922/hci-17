#ifndef HCI_TIMER_H_
#define HCI_TIMER_H_

#include <chrono>
#include <string>
#include <unordered_map>

#include "parallel.h"

class Time {
 private:
  static std::chrono::high_resolution_clock::time_point& get_time(const std::string& event) {
    using namespace std::chrono;
    static std::unordered_map<std::string, high_resolution_clock::time_point> times;
    return times[event];
  }

 public:
  static void init() {
    using namespace std::chrono;

    if (Parallel::get_id() != 0) return;
    Time::get_time("init") = high_resolution_clock::now();
  }
  static void start(const std::string& event) {
    using namespace std::chrono;

    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    const auto now = high_resolution_clock::now();
    const auto init_time = Time::get_time("init");
    Time::get_time(event) = now;
    const duration<double> since_init = duration_cast<duration<double>>(now - init_time);

    printf("[%.3f/0.000] BEG %s.\n", since_init.count(), event.c_str());
  };

  static void end(const std::string& event) {
    using namespace std::chrono;

    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    const auto now = high_resolution_clock::now();
    const auto start_time = Time::get_time(event);
    const auto init_time = Time::get_time("init");
    const duration<double> since_start = duration_cast<duration<double>>(now - start_time);
    const duration<double> since_init = duration_cast<duration<double>>(now - init_time);

    printf("[%.3f/%.3f] END %s.\n", since_init.count(), since_start.count(), event.c_str());
  }
};

#endif