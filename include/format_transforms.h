#pragma once
#include <Eigen/Geometry>

void format_transforms(const std::vector<Eigen::Affine3d>& T, std::string format, std::vector<Eigen::Affine3d>& A);