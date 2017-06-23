#include "helper_strings.h"

void HelperStrings::setup_ab() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    ab[dets[i].up.encode()].first.push_back(i);
    ab[dets[i].dn.encode()].second.push_back(i);
  }
  shrink(ab);
}

void HelperStrings::setup_ab_m1() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    const auto& up_elecs = dets[i].up.get_elec_orbs();
    SpinDet det_up(dets[i].up);
    for (std::size_t j = 0; j < up_elecs.size(); j++) {
      det_up.set_orb(up_elecs[j], false);
      ab_m1[det_up.encode()].first.push_back(i);
      det_up.set_orb(up_elecs[j], true);
    }

    const auto& dn_elecs = dets[i].dn.get_elec_orbs();
    SpinDet det_dn(dets[i].dn);
    for (std::size_t j = 0; j < dn_elecs.size(); j++) {
      det_dn.set_orb(dn_elecs[j], false);
      ab_m1[det_dn.encode()].second.push_back(i);
      det_dn.set_orb(dn_elecs[j], true);
    }
  }
  shrink(ab_m1);
}

void HelperStrings::shrink(
    std::unordered_map<Orbitals, std::pair<UnsignedInts, UnsignedInts>, boost::hash<Orbitals>>&
        helper_strings) {
  std::list<Orbitals> remove_list;
  for (auto& kv : helper_strings) {
    auto& value = kv.second;
    if (value.first.size() <= 1 && value.second.size() <= 1) remove_list.push_back(kv.first);
  }
  for (const auto& det : remove_list) helper_strings.erase(det);
}

UnsignedInts HelperStrings::find_potential_connections(const std::size_t i) {
  UnsignedInts connections;
  const Det& det = dets[i];
  const auto& up_elecs = det.up.get_elec_orbs();
  const auto& dn_elecs = det.dn.get_elec_orbs();

  // Add self.
  connections.push_back(i);
  connected[i] = true;

  // Two up/dn excitations.
  if (ab.find(det.dn.encode()) != ab.end()) {
    for (const std::size_t det_id : ab.find(det.dn.encode())->second.second) {
      if (!connected[det_id]) {
        connected[det_id] = true;
        connections.push_back(det_id);
      }
    }
  }
  if (ab.find(det.up.encode()) != ab.end()) {
    for (const std::size_t det_id : ab.find(det.up.encode())->second.first) {
      if (!connected[det_id]) {
        connected[det_id] = true;
        connections.push_back(det_id);
      }
    }
  }

  // One up one dn excitation.
  SpinDet det_up(det.up);
  SpinDet det_dn(det.dn);

  for (std::size_t i = 0; i < up_elecs.size(); i++) {
    det_up.set_orb(up_elecs[i], false);
    const auto& key_up = det_up.encode();
    const auto& kv_up = ab_m1.find(key_up);
    if (kv_up != ab_m1.end()) {
      for (const std::size_t det_id : kv_up->second.first) one_up[det_id] = true;
      for (std::size_t j = 0; j < dn_elecs.size(); j++) {
        det_dn.set_orb(dn_elecs[j], false);
        const auto& key_dn = det_dn.encode();
        if (ab_m1.find(key_dn) != ab_m1.end()) {
          for (const std::size_t det_id : ab_m1.find(key_dn)->second.second) {
            if (one_up[det_id] && !connected[det_id]) {
              connected[det_id] = true;
              connections.push_back(det_id);
            }
          }
        }
        det_dn.set_orb(dn_elecs[j], true);
      }
      for (const std::size_t det_id : kv_up->second.first) one_up[det_id] = false;
    }
    det_up.set_orb(up_elecs[i], true);
  }

  // Reset connected and return.
  for (const std::size_t det_id : connections) connected[det_id] = false;
  return connections;
}