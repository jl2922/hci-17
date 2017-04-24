#ifndef HCI_PARALLEL_H_
#define HCI_PARALLEL_H_

#include <boost/mpi.hpp>

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

  static std::string get_proc_name() { return Parallel::get_instance().env.processor_name(); }

  static void barrier() {
    fflush(stdout);
    Parallel::get_instance().world.barrier();
  }
};

#endif