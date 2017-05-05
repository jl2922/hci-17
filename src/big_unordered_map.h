#ifndef BIG_UNORDERED_MAP_H_
#define BIG_UNORDERED_MAP_H_

#ifndef SERIAL
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/utility.hpp>
#include <limits>
#include <vector>
#endif
#include <cstddef>
#include <unordered_map>
#include <utility>

#ifdef SERIAL
// Wrap normal hash map into the BigUnorderedMap interface.
template <class K, class V, class H>
class BigUnorderedMap {
 public:
  BigUnorderedMap(const std::pair<K, V>& skeleton) {}

  void reserve(unsigned long long n_buckets) { local_map.reserve(n_buckets); }

  unsigned long long bucket_count() const { return local_map.bucket_count(); }

  const std::unordered_map<K, V, H>& get_local_map() { return local_map; }

  void async_inc(const K& k, const V& v) { local_map[k] += v; }

  std::size_t get_target(const K&) { return 0; }

  void complete_async_incs(){};

  long long unsigned int size() const { return local_map.size(); }

 protected:
  std::unordered_map<K, V, H> local_map;
};
#else
// Distributed version for production run.
template <class K, class V, class H>
class BigUnorderedMap {
 public:
  BigUnorderedMap(const std::pair<K, V>& skeleton);

  void reserve(unsigned long long n_buckets);

  unsigned long long bucket_count() const;

  const std::unordered_map<K, V, H>& get_local_map() { return local_map; }

  void async_inc(const K&, const V&);

  std::size_t get_target(const K&);

  void complete_async_incs();

  long long unsigned int size() const;

 protected:
  // Basic MPI info.
  std::size_t proc_id;
  std::size_t n_procs;
  boost::mpi::communicator world;

  // Hash function.
  H hasher;

  // Local hash map.
  std::unordered_map<K, V, H> local_map;

  // Buffers.
  std::vector<std::pair<K, V> > buf_send;
  std::vector<boost::mpi::content> buf_send_content;
  std::pair<K, V> buf_recv;
  boost::mpi::content buf_recv_content;
  std::list<boost::mpi::request> reqs;
  std::size_t buf_send_cnt;
  std::size_t buf_size;

  // Counters.
  std::vector<unsigned long long> send_cnts;
  std::vector<unsigned long long> recv_cnts;
  std::vector<unsigned long long> recv_totals;
  std::vector<unsigned long long> recv_trunk_totals;

  // MPI tags.
  enum { TAG_NODE_INFO, TAG_KV, TAG_TRUNK_FINISH, TAG_FINISH };

  // Proc infos.
  std::size_t total_proc_buckets;
  std::size_t local_proc_buckets;
  std::vector<std::size_t> proc_map;

  void set_skeleton(const std::pair<K, V>&);

  void reset_cnts();

  void set_proc_buckets();

  void complete_async_inc_trunk();
};

template <class K, class V, class H>
BigUnorderedMap<K, V, H>::BigUnorderedMap(const std::pair<K, V>& skeleton) {
  const std::size_t TOTAL_BUF_SIZE = 2000;
  const std::size_t MIN_BUF_SIZE = 100;
  proc_id = world.rank();
  n_procs = world.size();
  hasher = H();
  buf_size = std::max(TOTAL_BUF_SIZE / n_procs, MIN_BUF_SIZE);
  buf_send.resize(buf_size);
  buf_send_content.resize(buf_size);
  reset_cnts();
  set_skeleton(skeleton);
  set_proc_buckets();
}

template <class K, class V, class H>
unsigned long long BigUnorderedMap<K, V, H>::bucket_count() const {
  unsigned long long local_bucket_count = static_cast<unsigned long long>(local_map.bucket_count());
  unsigned long long total_bucket_count = 0;
  boost::mpi::all_reduce(
      world, local_bucket_count, total_bucket_count, std::plus<unsigned long long>());
  return total_bucket_count;
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::set_skeleton(const std::pair<K, V>& skeleton) {
  for (std::size_t i = 0; i < buf_size; i++) {
    buf_send[i] = skeleton;
    buf_send_content[i] = boost::mpi::get_content(buf_send[i]);
  }
  buf_recv = skeleton;
  buf_recv_content = boost::mpi::get_content(buf_recv);
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::reset_cnts() {
  buf_send_cnt = 0;
  send_cnts.assign(n_procs, 0);
  recv_cnts.assign(n_procs, 0);
  recv_totals.assign(n_procs, std::numeric_limits<unsigned long long>::max());
  recv_trunk_totals.assign(n_procs, std::numeric_limits<unsigned long long>::max());
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::set_proc_buckets() {
  // Setup proc buckets according the memory size stored in nodes.json.
  // Proc buckets controls the storage load balance.
  const std::size_t NODE_SIZE_DEFAULT = 8;
  const std::size_t TOTAL_INTERNAL_UNITS = 10000;
  boost::property_tree::ptree nodes;
  boost::property_tree::read_json("nodes.json", nodes);
  std::vector<std::size_t> node_sizes(n_procs, 0);  // Memory size per node.
  std::vector<std::size_t> proc_sizes(n_procs, 0);  // Memory size (in internal units) per process.
  std::vector<std::string> node_names(n_procs);
  std::unordered_map<std::string, std::size_t> node_procs;

  // Send local node name to other processes.
  std::string local_node_name = boost::mpi::environment::processor_name();
  std::list<boost::mpi::request> reqs;
  for (std::size_t i = 0; i < n_procs; i++) {
    reqs.push_front(world.isend(i, TAG_NODE_INFO, local_node_name));
  }

  // Receive node buckets info from remote nodes.
  std::size_t total_size = 0;
  std::size_t recv_cnt = 0;
  while (recv_cnt < n_procs) {
    const auto& status = world.probe(boost::mpi::any_source, TAG_NODE_INFO);
    const int source = status.source();
    std::string node_name;
    world.recv(source, TAG_NODE_INFO, node_name);
    node_names[source] = node_name;
    node_sizes[source] = nodes.get<std::size_t>(
        boost::property_tree::ptree::path_type(node_name, '/'), NODE_SIZE_DEFAULT);
    total_size += node_sizes[source];
    node_procs[node_name]++;
    recv_cnt++;
  }
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  // Obtain process sizes.
  total_proc_buckets = 0;
  proc_map.clear();
  proc_map.reserve(TOTAL_INTERNAL_UNITS);
  for (std::size_t i = 0; i < n_procs; i++) {
    proc_sizes[i] = node_sizes[i] * TOTAL_INTERNAL_UNITS / total_size / node_procs[node_names[i]];
    total_proc_buckets += proc_sizes[i];
    for (std::size_t j = 0; j < proc_sizes[i]; j++) proc_map.push_back(i);
  }
  local_proc_buckets = proc_sizes[proc_id];

  // Print balance info.
  printf(
      "Proc #%lu (on %s) stores: %.3f %%\n",
      proc_id,
      local_node_name.c_str(),
      local_proc_buckets * 100.0 / total_proc_buckets);
  world.barrier();
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::reserve(const unsigned long long n_buckets) {
  std::size_t local_buckets =
      static_cast<std::size_t>(n_buckets * local_proc_buckets / total_proc_buckets + 1);
  local_map.reserve(local_buckets);
  world.barrier();
}

template <class K, class V, class H>
std::size_t BigUnorderedMap<K, V, H>::get_target(const K& key) {
  return proc_map[hasher(key) % total_proc_buckets];
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::async_inc(const K& key, const V& value) {
  const std::size_t target = get_target(key);

  // Process locally.
  if (target == proc_id) {
    local_map[key] += value;
    send_cnts[proc_id]++;
    return;
  }

  // Send to target asynchronously.
  buf_send[buf_send_cnt].first = key;
  buf_send[buf_send_cnt].second = value;
  reqs.push_front(world.isend(target, TAG_KV, buf_send_content[buf_send_cnt]));
  send_cnts[target]++;
  if (send_cnts[target] == std::numeric_limits<unsigned long long>::max()) {
    printf("Maximum sends reached. Call complete async first.\n");
    exit(1);
  }
  buf_send_cnt++;

  if (buf_send_cnt == buf_size) {
    // When buffer is full, wait and process a trunk.
    // A trunk contains 'buf_size' send-requests from each process.
    for (std::size_t i = 0; i < n_procs; i++) {
      if (i == proc_id) continue;
      reqs.push_front(world.isend(i, TAG_TRUNK_FINISH, send_cnts[i]));
    }
    complete_async_inc_trunk();
  }
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::complete_async_inc_trunk() {
  std::size_t n_active_procs = 0;
  std::size_t trunk_finish_cnt = 0;
  for (std::size_t i = 0; i < n_procs; i++) {
    if (recv_cnts[i] < recv_totals[i] && i != proc_id) n_active_procs++;
  }
  while (trunk_finish_cnt < n_active_procs) {
    const auto& status = world.probe();
    const std::size_t source = status.source();
    const std::size_t tag = status.tag();
    switch (tag) {
      case TAG_KV: {
        world.recv(source, tag, buf_recv_content);
        const K& key = buf_recv.first;
        const V& value = buf_recv.second;
        local_map[key] += value;
        recv_cnts[source]++;
        if (recv_cnts[source] == recv_trunk_totals[source]) {
          trunk_finish_cnt++;
        }
        break;
      }
      case TAG_TRUNK_FINISH: {
        world.recv(source, TAG_TRUNK_FINISH, recv_trunk_totals[source]);
        if (recv_trunk_totals[source] == recv_cnts[source]) {
          trunk_finish_cnt++;
        }
        break;
      }
      case TAG_FINISH: {
        world.recv(source, TAG_FINISH, recv_totals[source]);
        if (recv_totals[source] == recv_cnts[source]) trunk_finish_cnt++;
        break;
      }
    }
  }
  wait_all(reqs.begin(), reqs.end());
  reqs.clear();
  buf_send_cnt = 0;
  recv_trunk_totals.assign(n_procs, std::numeric_limits<unsigned long long>::max());
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::complete_async_incs() {
  std::size_t n_active_procs = 0;
  for (std::size_t i = 0; i < n_procs; i++) {
    if (i == proc_id) continue;
    reqs.push_front(world.isend(i, TAG_FINISH, send_cnts[i]));
    if (recv_cnts[i] < recv_totals[i]) n_active_procs++;
  }

  while (n_active_procs > 0) {
    const auto& status = world.probe();
    const int source = status.source();
    const int tag = status.tag();
    switch (tag) {
      case TAG_KV: {
        world.recv(source, tag, buf_recv_content);
        const K& key = buf_recv.first;
        const V& value = buf_recv.second;
        local_map[key] += value;
        recv_cnts[source]++;
        if (recv_cnts[source] == recv_totals[source]) n_active_procs--;
        break;
      }
      case TAG_TRUNK_FINISH: {
        world.recv(source, TAG_TRUNK_FINISH, recv_trunk_totals[source]);
        break;
      }
      case TAG_FINISH: {
        world.recv(source, TAG_FINISH, recv_totals[source]);
        if (recv_totals[source] == recv_cnts[source]) n_active_procs--;
        break;
      }
    }
  }
  wait_all(reqs.begin(), reqs.end());
  reqs.clear();
  reset_cnts();
  world.barrier();
}

template <class K, class V, class H>
unsigned long long BigUnorderedMap<K, V, H>::size() const {
  unsigned long long local_size = static_cast<unsigned long long>(local_map.size());
  unsigned long long total_size = 0;
  boost::mpi::all_reduce(world, local_size, total_size, std::plus<unsigned long long>());
  return total_size;
}
#endif

#endif