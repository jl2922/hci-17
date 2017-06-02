#ifndef TIMER_H_
#define TIMER_H_

#include "std.h"

class Timer {
 public:
  static void now() {
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << std::ctime(&t);
  }

  static void start(const std::string& event) {
    Timer& timer = Timer::get_instance();
    timer.timepoints.push_back(
        std::pair<std::string, std::chrono::high_resolution_clock::time_point>(
            event, std::chrono::high_resolution_clock::now()));
    timer.next_indices.push_back(0);
    printf("\n%sBEG ", timer.get_current_indices().c_str());
    printf("%s", timer.get_current_event_path().c_str());
    printf(" [%.3f]\n", get_duration(timer.init_time));
  }

  static void end() {
    Timer& timer = Timer::get_instance();
    printf("%sEND ", timer.get_current_indices().c_str());
    printf("%s", timer.get_current_event_path().c_str());
    printf(
        " [%.3f/%.3f]\n",
        get_duration(timer.init_time),
        get_duration(timer.timepoints.back().second));
    timer.timepoints.pop_back();
    timer.next_indices.pop_back();
    timer.next_indices.back()++;
  }

 private:
  std::vector<std::pair<std::string, std::chrono::high_resolution_clock::time_point>> timepoints;
  std::vector<std::size_t> next_indices;
  std::chrono::high_resolution_clock::time_point init_time;

  Timer() {
    init_time = std::chrono::high_resolution_clock::now();
    next_indices.push_back(0);
  }

  std::string get_current_event_path() {
    std::stringstream ss;
    for (std::size_t i = 0; i < timepoints.size(); i++) {
      ss << timepoints[i].first;
      if (i != timepoints.size() - 1) {
        ss << " >> ";
      }
    }
    return ss.str();
  }

  std::string get_current_indices() {
    std::stringstream ss;
    for (std::size_t i = 1; i < timepoints.size(); i++) {
      ss << next_indices[i] << '.';
    }
    return ss.str();
  }

  static double get_duration(std::chrono::high_resolution_clock::time_point start) {
    const auto now = std::chrono::high_resolution_clock::now();
    return (std::chrono::duration_cast<std::chrono::duration<double>>(now - start)).count();
  }

  static Timer& get_instance() {
    static Timer timer;
    return timer;
  }
};

#endif