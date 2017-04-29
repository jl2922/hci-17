#ifndef HCI_PARALLEL_H_
#define HCI_PARALLEL_H_

#include <boost/mpi.hpp>
#include "std.h"

class Parallel {
 private:
  int id;
  int n;
  boost::mpi::environment env;
  boost::mpi::communicator world;

  Parallel() {
    this->id = this->world.rank();
    this->n = this->world.size();
  }

  // Singleton pattern.
  static Parallel& get_instance() {
    static Parallel instance;
    return instance;
  }

 public:
  static int get_id() { return Parallel::get_instance().id; }

  static int get_n() { return Parallel::get_instance().n; }

  static std::string get_host() { return Parallel::get_instance().env.processor_name(); }

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

#endif