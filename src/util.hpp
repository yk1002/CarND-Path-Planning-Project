#pragma once

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

// Debug function to have the compiler show the type of T
template<typename T> void showtype(T&& parameter);  // purposefully not defined

// Math utility functions
inline constexpr double pi() { return M_PI; }
inline double deg2rad(double x) { return x * pi() / 180; }
inline double rad2deg(double x) { return x * 180 / pi(); }
inline double distance(double x1, double y1, double x2, double y2) {
  const auto xd = x1 - x2;
  const auto yd = y1 - y2;
  return sqrt(xd * xd + yd * yd);
}

// Frame conversion
Eigen::Matrix3d getLocal2GlobalConvMatrix(double px, double py, double psi);
Eigen::Matrix3d getGlobal2LocalConvMatrix(double px, double py, double psi);

template <typename T1, typename T2>
void convertFrame(const Eigen::Matrix3d& conversionMatrix,
                  const T1& ptsxFrom, const T1& ptsyFrom, T2& ptsxTo, T2& ptsyTo) {
  const int N_pts = ptsxFrom.size();
  assert(N_pts == ptsxTo.size());
  assert(N_pts == ptsyTo.size());
  for (int i = 0; i < N_pts; ++i) {
    const Eigen::Vector3d ptsFrom(ptsxFrom[i], ptsyFrom[i], 1.0);
    const auto ptsTo = conversionMatrix * ptsFrom;
    ptsxTo[i] = ptsTo[0];
    ptsyTo[i] = ptsTo[1];
  }
}

// For printing vector to ostream
template <typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& container) {
  auto i = begin(container);
  const auto iEnd = end(container);
  os << "{";
  bool first = true;
  while (i != iEnd) {
    os << (first ? "" : ", " ) << *i;
    first = false;;
    ++i;
  }
  os << "}";
  return os;
}

// For converting Frenet coordinates to Cartesian coordinates and vice versa
class MapUtil {
 public:
  const double max_s;
  
  /**
     map_file: map file path
     max_s: s value before wrapping around the track back to 0
  */
  MapUtil(std::string map_file, double max_s);

  std::pair<double, double> waypoint(int i) const { return std::make_pair(maps_x[i], maps_y[i]); }
  int closestWaypointIndex(double x, double y) const;
  int nextWaypointIndex(double x, double y, double theta) const;

  // Transform from Cartesian (x, y) coordinates to Frenet (s, d) coordinates
  std::pair<double, double> getFrenet(double x, double y, double theta) const;

  // Transform from Frenet (s, d) coordinates to Cartesian (x, y) coordinates
  std::pair<double, double> getXY(double s, double d) const;

 private:
  std::vector<double> maps_x;
  std::vector<double> maps_y;
  std::vector<double> maps_s;
  std::vector<double> maps_dx;
  std::vector<double> maps_dy;
  int nWaypoints;
};
