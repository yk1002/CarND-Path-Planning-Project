#include "util.hpp"
#include <fstream>

using namespace std;

/* local to global conversion:
   | Xg |    | cos(psi) -sin(psi) px | | Xc |
   | Yg | =  | sin(psi)  cos(psi) py | | Yc |
   | 1  |    |        0         0  1 | |  1 |
*/
Eigen::Matrix3d getLocal2GlobalConvMatrix(double px, double py, double psi) {
  const auto c = cos(psi);
  const auto s = sin(psi);
  Eigen::Matrix3d m;
  m <<
      c,    -s,  px,
      s,     c,  py,
      0.0, 0.0, 1.0;
  return m;
}

/* global to local conversion (inverse of the matrix above):
   | Xc |    |  cos(psi)  sin(psi) -sin(psi)*py-cos(psi)*px | | Xg |
   | Yc | =  | -sin(psi)  cos(psi)  sin(psi)*px-cos(psi)*py | | Yg |
   | 1  |    |         0         0                        1 | |  1 |
*/
Eigen::Matrix3d getGlobal2LocalConvMatrix(double px, double py, double psi) {
  const auto c = cos(psi);
  const auto s = sin(psi);
  Eigen::Matrix3d m;
  m << 
      c,     s, -s * py - c * px,
      -s,    c,  s * px - c * py,
      0.0, 0.0,              1.0;
  return m;
}

MapUtil::MapUtil(string map_file, double max_s) : max_s(max_s) {
  // Load up map values for waypoint's x, y, s and d normalized normal vectors
  ifstream ifs(map_file.c_str(), ifstream::in);

  double x, y, s, dx, dy;
  while (ifs >> x >> y >> s >> dx >> dy) {
    maps_x.push_back(x);
    maps_y.push_back(y);
    maps_s.push_back(s);
    maps_dx.push_back(dx);
    maps_dy.push_back(dy);
  }

  nWaypoints = maps_x.size();
}

int MapUtil::closestWaypointIndex(double x, double y) const {
  double closestDist = 1e60;
  int closestIndex = 0;
  
  for (int i = 0; i < nWaypoints; ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestDist) {
      closestDist = dist;
      closestIndex = i;
    }
  }

  return closestIndex;
}

int MapUtil::nextWaypointIndex(double x, double y, double theta) const {
  auto nextIndex = closestWaypointIndex(x, y);
  double map_x = maps_x[nextIndex];
  double map_y = maps_y[nextIndex];
  double heading = atan2(map_y - y, map_x - x);
  double angle = abs(theta - heading);
  if (angle > pi()/4) ++nextIndex;
  return nextIndex;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
pair<double, double> MapUtil::getFrenet(double x, double y, double theta) const {
  int next_wp = nextWaypointIndex(x, y, theta);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = nWaypoints-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return make_pair(frenet_s, frenet_d);
}

// Transform from Frenet s,d coordinates to Cartesian x,y
pair<double, double> MapUtil::getXY(double s, double d) const {
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%nWaypoints;

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return make_pair(x, y);
}
