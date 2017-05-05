#ifndef HCI_TIME_H_
#define HCI_TIME_H_

#include "../std.h"

#include "../parallel.h"
#include "index.h"

class Time {
 public:
  static void init() {
    if (Parallel::get_id() != 0) return;
    Time::get_timer("INIT") = std::chrono::high_resolution_clock::now();
  }

  static void start(const std::string& event) {
    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    Time::get_timer(event) = std::chrono::high_resolution_clock::now();
    auto& index = Time::get_index();
    index.start();
    printf(
        "\n%sBEGIN %s [%.3f/0.000]\n",
        index.to_string().c_str(),
        event.c_str(),
        get_duration("INIT"));
  };

  static void end(const std::string& event) {
    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    auto& index = Time::get_index();
    printf(
        "%sEND %s [%.3f/%.3f] \n",
        index.to_string().c_str(),
        event.c_str(),
        get_duration("INIT"),
        get_duration(event));
    index.end();
  }

  static void checkpoint(const std::string& event) {
    printf(
        "CHECKPOINT %s [%.3f/%.3f] \n", event.c_str(), get_duration("INIT"), get_duration(event));
  }

 private:
  static Index& get_index() {
    static Index index;
    return index;
  };

  static std::chrono::high_resolution_clock::time_point& get_timer(const std::string& event) {
    static std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> times;
    return times[event];
  }

  static double get_duration(const std::string& event) {
    using namespace std::chrono;
    const auto now = high_resolution_clock::now();
    const auto start_time = Time::get_timer(event);
    return (duration_cast<duration<double>>(now - start_time)).count();
  }
};

#endif