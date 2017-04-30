#ifndef HCI_PARALLEL_H_
#define HCI_PARALLEL_H_

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "std.h"

#ifndef SERIAL
class Parallel {
 private:
  int id;
  int n;
  boost::mpi::environment* env;  // For MPI 1.1.
  boost::mpi::communicator world;

  Parallel() {
    id = this->world.rank();
    n = this->world.size();
  }

  // Singleton pattern.
  static Parallel& get_instance() {
    static Parallel instance;
    return instance;
  }

 public:
  static void init(boost::mpi::environment& env) { Parallel::get_instance().env = &env; }

  static int get_id() { return Parallel::get_instance().id; }

  static int get_n() { return Parallel::get_instance().n; }

  static std::string get_host() { return Parallel::get_instance().env->processor_name(); }

  static void barrier() {
    fflush(stdout);
    Parallel::get_instance().world.barrier();
  }

  template <class T>
  static void all_reduce(T& t) {
    T t_local = t;
    boost::mpi::all_reduce(Parallel::get_instance().world, t_local, t, std::plus<T>());
  }
};
#else
// Non-MPI stub for debugging and memory profiling.
class Parallel {
 private:
  Parallel() {}

  // Singleton pattern.
  static Parallel& get_instance() {
    static Parallel instance;
    return instance;
  }

 public:
  static int get_id() { return 0; }

  static int get_n() { return 1; }

  static std::string get_host() { return "localhost"; }

  static void barrier() { return; }

  template <class T>
  static void all_reduce(T& t) {
    return;
  }
};
#endif

#endif