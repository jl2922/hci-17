#ifndef BIG_UNORDERED_MAP_H_
#define BIG_UNORDERED_MAP_H_

#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

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

#endif