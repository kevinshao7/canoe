#ifndef SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
#define SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

namespace CubedSphereUtility {

int encode_offset(int pid, int ox2, int ox3);

int decode_offset(int *pid, int *ox2, int *ox3, int offset);

int find_panel_id(LogicalLocation const &loc) {
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  return lv2_lx2 + lv2_lx3 * 2 + 1;
}

void find_local_position(int *local_lx2, int *local_lx3, int *bound_lim,
                         LogicalLocation const &loc) {
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);

  *local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  *local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
  *bound_lim = (1 << (loc.level - 2)) - 1;
}

int find_target_offset(int ox2, int ox3, LogicalLocation const &loc) {
  int pid = find_panel_id(loc);
  int offset = encode_offset(pid, ox2, ox3);
  return target_offset_[offset];
}

int find_source_offset(int tox2, int tox3, LogicalLocation const &loc) {
  int pid = find_panel_id(loc);
  int offset = encode_offset(pid, tox2, tox3);
  return source_offset_[offset];
}

Real generate_mesh_x2(Real x, LogicalLocation const &loc);

Real generate_mesh_x3(Real x, LogicalLocation const &loc);

void PackData(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
              int ei, int sj, int ej, int sk, int ek, int &offset, int ox1,
              int ox2, int ox3, LogicalLocation const &loc);

// Helper functions adapted from Paul
void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3);
void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V);
void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat);
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP);

template <typename A>
void CovariantToContravariant(A a, Real cth) {
  Real v = a[IVY];
  Real w = a[IVZ];
  Real sth2 = 1. - cth * cth;

  a[IVY] = v / sth2 - w * cth / sth2;
  a[IVZ] = -v * cth / sth2 + w / sth2;
}

template <typename A>
void ContravariantToCovariant(A a, Real cth) {
  Real v = a[IVY];
  Real w = a[IVZ];
  a[IVY] = v + w * cth;
  a[IVZ] = w + v * cth;
}

}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
