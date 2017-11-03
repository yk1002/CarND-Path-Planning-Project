#include "PathPlanner.hpp"

#include <iomanip>
#include <iostream>
#include "spline.h"

bool verbose = true;
//bool verbose = false;
#define V(v) #v": "<< v

using namespace std;

PathPlanner::PathPlanner(const MapUtil& mapUtil) :
    mapUtil(mapUtil) {
}

  // // resuse previous points
  // copy(begin(previous_path_x), end(previous_path_x), back_inserter(path_x));
  // copy(begin(previous_path_y), end(previous_path_y), back_inserter(path_y));

/*
  TODO How to make a decision on lane change? (keep lane or shift to left or right)
  TODO How to decide where to aim for?
  TODO How to avoid collision (side, front, back)?
  TODO How to use FSM?
 */
pair<vector<double>, vector<double>>
PathPlanner::calculatePath(double car_x, double car_y, double car_s, double car_d, double car_yaw, double car_speed,
                           const vector<double>& previous_path_x, const vector<double>& previous_path_y, 
                           double end_path_s, double end_path_d, const vector<vector<double>>& sensor_fusion) {
  if (verbose) {
    static int seq = 0;
    cout << setprecision(5);
    cout << "===== seq: " << seq++ << "=====" << '\n'
         << V(car_x) << '\n'
         << V(car_y) << '\n'
         << V(car_s) << '\n'
         << V(car_d) << '\n'
         << V(car_yaw) << '\n'
         << V(car_speed) << '\n'
         << V(previous_path_x) << '\n'
         << V(previous_path_y) << '\n'
         << V(end_path_s) << '\n'
         << V(end_path_d) << endl;

    //    if (seq == 3) exit(1);
  }
  
  const int prev_size = previous_path_x.size();

  // Add spline anchor points in global coordinates.
  vector<double> splineGx, splineGy;
  double splineX0, splineY0, splineX1, splineY1, splineStartYaw; 
  if (prev_size >= 2) { // first two points are picked from previous path
    splineX0 = previous_path_x[prev_size - 2];
    splineY0 = previous_path_y[prev_size - 2];
    splineX1 = previous_path_x[prev_size - 1];
    splineY1 = previous_path_y[prev_size - 1];
  }
  else {               // ...or current position if previous path is absent
    splineX0 = car_x - cos(car_yaw);
    splineY0 = car_y - sin(car_yaw);
    splineX1 = car_x;
    splineY1 = car_y;
  }
  splineStartYaw = atan2(splineY1 - splineY0, splineX1 - splineX0);
  splineGx.push_back(splineX0);
  splineGy.push_back(splineY0);
  splineGx.push_back(splineX1);
  splineGy.push_back(splineY1);

  int lane = 1;
  const auto target_d = 2.0 + lane * 4;
  for (auto dist: { 30, 60, 90 }) { // in meters; 50 [mph] = 20 [m/s]
    const auto p = mapUtil.getXY(car_s + dist, target_d);
    splineGx.push_back(p.first);
    splineGy.push_back(p.second);
  }
  
  // Convert spline anchor points to local coordinates
  vector<double> splineLx(splineGx.size()), splineLy(splineGx.size());
  const auto g2l = getGlobal2LocalConvMatrix(splineX1, splineY1, splineStartYaw);
  convertFrame(g2l, splineGx, splineGy, splineLx, splineLy);

  // check splineLx values are monotonically increasing
  for (int i = 1; i < splineLx.size(); ++i) {
    if (splineLx[i - 1] >= splineLx[i]) {
      cout << "\n============== ERROR! ==============" << '\n'
           << "splineLx[" << i-1 << "]=" << splineLx[i - 1] << '\n'
           << "splineLx[" << i << "]=" << splineLx[i] << '\n'
           << V(splineX0) << '\n'
           << V(splineY0) << '\n'
           << V(splineStartYaw) << '\n'
           << V(splineGx) << '\n'
           << V(splineGy) << '\n'
           << V(splineLx) << '\n'
           << V(splineLy) << '\n'
           << "\n============== ERROR! ==============" << endl;
    }
  }

  // Fit spline curve over anchor points
  tk::spline sp;
  sp.set_points(splineLx, splineLy);

  // Now, pick point on the spline curve in such a way that
  // picked points are spaced "legally."
  const auto timeIntervalS = 0.02;  // seconds
  const auto maxSpeedMPH = 50.0; // miles per hour
  const auto speedMPS = car_speed * 1.60934 * 1000 / 3600; // meters per second
  const auto maxSpeedMPS = maxSpeedMPH * 1.60934 * 1000 / 3600; // meters per second
  const auto maxIntervalM = maxSpeedMPS * timeIntervalS; // max distance allowed between points in meters

  // Now, sample points on the spline in such a way that
  // 1. Any two adjecent points are less than maxIntervalM apart.
  // 2. Acceleration does not exceed (TBD)
  const int N = 50;
  const int nPoints = N - prev_size; // 0.02 [s] * 50 = 1.0 [s]
  vector<double> lx(nPoints), ly(nPoints);
  double prevX = 0;
  double prevY = sp(prevX);
  for (int i = 0; i < nPoints; ++i) {
    auto x = prevX + maxIntervalM;
    auto y = sp(x);
    const auto slope = distance(prevX, prevY, x, y);
    const auto factor = maxIntervalM / slope * 0.9;
    x = prevX + maxIntervalM * factor;
    y = sp(x);
    lx[i] = x;
    ly[i] = y;
    prevX = x;
    prevY = y;
  }

  // Convert sample points in local coordinates to global coordinates.
  vector<double> gx(lx.size()), gy(ly.size());
  const auto l2g = getLocal2GlobalConvMatrix(splineX1, splineY1, splineStartYaw);
  convertFrame(l2g, lx, ly, gx, gy);

  vector<double> ptx, pty;
  copy(begin(previous_path_x), end(previous_path_x), back_inserter(ptx));
  copy(begin(previous_path_y), end(previous_path_y), back_inserter(pty));
  copy(begin(gx), end(gx), back_inserter(ptx));
  copy(begin(gy), end(gy), back_inserter(pty));
  assert(ptx.size() == N);

  if (verbose) {
    cout << V(splineGx) << '\n'
         << V(splineGy) << '\n'
         << V(splineLx) << '\n'
         << V(splineLy) << '\n'
         << V(lx) << '\n'
         << V(ly) << '\n'
         << V(gx) << '\n'
         << V(gy) << '\n'
         << V(ptx) << '\n'
         << V(pty) << endl;
  }

  return make_pair(ptx, pty);
}
