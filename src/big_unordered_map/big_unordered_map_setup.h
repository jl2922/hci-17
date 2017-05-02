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
  const std::size_t NODE_SIZE_DEFAULT = 8;
  const std::size_t TOTAL_INTERNAL_UNITS = 10000;
  boost::property_tree::ptree nodes;
  boost::property_tree::read_json("nodes.json", nodes);
  std::vector<std::size_t> node_sizes(n_procs, 0);  // Memory size per node.
  std::vector<std::size_t> proc_sizes(n_procs, 0);  // Memory size (in internal unit) per process.
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
      static_cast<std::size_t>(n_buckets * 1.0 * local_proc_buckets / total_proc_buckets + 1);
  local_map.reserve(local_buckets);
  world.barrier();
}