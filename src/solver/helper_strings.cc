#include "helper_strings.h"

void HelperStrings::setup_ab() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    const auto& up_elecs = dets[i].up.get_elec_orbs();
    ab[up_elecs].first.push_back(i);

    const auto& dn_elecs = dets[i].dn.get_elec_orbs();
    ab[dn_elecs].second.push_back(i);
  }
  shrink(ab);
}

void HelperStrings::setup_ab_m1() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    const auto& up_elecs = dets[i].up.get_elec_orbs();
    SpinDet det_up(dets[i].up);
    for (std::size_t j = 0; j < up_elecs.size(); j++) {
      det_up.set_orb(up_elecs[j], false);
      ab_m1[det_up.get_elec_orbs()].first.push_back(i);
      det_up.set_orb(up_elecs[j], true);
    }

    const auto& dn_elecs = dets[i].dn.get_elec_orbs();
    SpinDet det_dn(dets[i].dn);
    for (std::size_t j = 0; j < dn_elecs.size(); j++) {
      det_dn.set_orb(dn_elecs[j], false);
      ab_m1[det_dn.get_elec_orbs()].second.push_back(i);
      det_dn.set_orb(dn_elecs[j], true);
    }
  }
  shrink(ab_m1);
}

void HelperStrings::shrink(
    std::unordered_map<Orbitals, std::pair<Ints, Ints>, boost::hash<Orbitals>>& helper_strings) {
  std::list<Orbitals> remove_list;
  for (auto& kv : helper_strings) {
    auto& value = kv.second;
    if (value.first.size() <= 1 && value.second.size() <= 1) {
      remove_list.push_back(kv.first);
      continue;
    }
  }
  for (const auto& det : remove_list) helper_strings.erase(det);
}

Ints HelperStrings::find_potential_connections(const int i) {
  Ints connections;
  const Det& det = dets[i];
  const auto& up_elecs = det.up.get_elec_orbs();
  const auto& dn_elecs = det.dn.get_elec_orbs();

  // Add self.
  connections.push_back(i);
  connected[i] = true;

  // Two up/dn excitations.
  if (ab.find(dn_elecs) != ab.end()) {
    for (const int orb_id : ab.find(dn_elecs)->second.second) {
      if (!connected[orb_id]) {
        connected[orb_id] = true;
        connections.push_back(orb_id);
      }
    }
  }
  if (ab.find(up_elecs) != ab.end()) {
    for (const int orb_id : ab.find(up_elecs)->second.first) {
      if (!connected[orb_id]) {
        connected[orb_id] = true;
        connections.push_back(orb_id);
      }
    }
  }

  // One up one dn excitation.
  SpinDet det_up(det.up);
  SpinDet det_dn(det.dn);

  for (std::size_t i = 0; i < up_elecs.size(); i++) {
    det_up.set_orb(up_elecs[i], false);
    const auto& key_up = det_up.get_elec_orbs();
    if (ab_m1.find(key_up) != ab_m1.end()) {
      for (const int orb_id : ab_m1.find(key_up)->second.first) one_up[orb_id] = true;
      for (std::size_t j = 0; j < dn_elecs.size(); j++) {
        det_dn.set_orb(dn_elecs[j], false);
        const auto& key_dn = det_dn.get_elec_orbs();
        if (ab_m1.find(key_dn) != ab_m1.end()) {
          for (const int orb_id : ab_m1.find(key_dn)->second.second) {
            if (one_up[orb_id] && !connected[orb_id]) {
              connected[orb_id] = true;
              connections.push_back(orb_id);
            }
          }
        }
        det_dn.set_orb(dn_elecs[j], true);
      }
      for (const int orb_id : ab_m1.find(key_up)->second.first) one_up[orb_id] = false;
    }
    det_up.set_orb(up_elecs[i], true);
  }

  // Reset connected and return.
  for (const int orb_id : connections) connected[orb_id] = false;  // Prepare for the next call.
  return connections;
}