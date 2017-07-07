#ifndef HCI_TIME_H_
#define HCI_TIME_H_

#include "../std.h"

#include "../parallel.h"
#include "index.h"

class Time {
 public:
  // To be called at the beginning of the program.
  static void init() {
    if (Parallel::get_id() != 0) return;
    Time::get_instance().timers["INIT"] = std::chrono::high_resolution_clock::now();
  }

  // To be called at the start of an event.
  static void start(const std::string& event) {
    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    Time::get_instance().timers[event] = std::chrono::high_resolution_clock::now();
    auto& index = Time::get_instance().index;
    index.start();
    printf(
        "\n%sBEGIN %s [%.3f/0.000]\n",
        index.to_string().c_str(),
        event.c_str(),
        get_duration("INIT"));
  };

  // To be called at the end of an event.
  static void end(const std::string& event) {
    Parallel::barrier();
    if (Parallel::get_id() != 0) return;
    auto& index = Time::get_instance().index;
    printf(
        "%sEND %s [%.3f/%.3f] \n",
        index.to_string().c_str(),
        event.c_str(),
        get_duration("INIT"),
        get_duration(event));
    index.end();
  }

  // To be called during an event after an important step.
  static void checkpoint(const std::string& event, const std::string& msg = "") {
    if (Parallel::get_id() != 0) return;
    printf(
        "CHECKPOINT %s: %s [%.3f/%.3f] \n",
        event.c_str(),
        msg.c_str(),
        get_duration("INIT"),
        get_duration(event));
  }

 private:
  Index index;
  std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> timers;

  Time() {}

  // Singleton pattern.
  static Time& get_instance() {
    static Time instance;
    return instance;
  }

  static double get_duration(const std::string& event) {
    using namespace std::chrono;
    const auto now = high_resolution_clock::now();
    const auto start_time = Time::get_instance().timers[event];
    return (duration_cast<duration<double>>(now - start_time)).count();
  }
};

#endif