#pragma once
#include <Eigen/Core>

bool write_anim_to_json(const Eigen::MatrixXd& P, const std::string file, const std::string format = "blender");