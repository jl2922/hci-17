#ifndef HCI_INDEX_H_
#define HCI_INDEX_H_

#include "../std.h"

class Index {
 public:
  Index() { index.push_back(0); }

  // Print out x.x.x.x.
  std::string to_string() {
    const int INDEX_IGNORE = 2;
    if (INDEX_IGNORE >= index.size() - 1) return "";
    std::stringstream ss;
    for (std::size_t i = INDEX_IGNORE; i < index.size() - 1; i++) ss << index[i] << '.';
    return ss.str() + " ";
  }

  // Add sub index under current index.
  void start() {
    index.back()++;
    index.push_back(0);
  }

  // Finish current sub index.
  void end() { index.pop_back(); }

 private:
  std::vector<int> index;
};

#endif