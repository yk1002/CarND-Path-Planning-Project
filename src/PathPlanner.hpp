#pragma once

#include <vector>
#include "util.hpp"

class PathPlanner {
 public:
  PathPlanner(const MapUtil& mapUtil);

  std::pair<std::vector<double>, std::vector<double>>
  calculatePath(double car_x, double car_y, double car_s, double car_d, double car_yaw, double car_speed,
                const std::vector<double>& previous_path_x, const std::vector<double>& previous_path_y, 
                double end_path_s, double end_path_d, const std::vector<std::vector<double>>& sensor_fusion);

 private:
  const MapUtil& mapUtil;
};
