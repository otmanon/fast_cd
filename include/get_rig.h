#pragma once
#include "rig.h"


Rig* get_rig(std::string rig_file, Eigen::MatrixXd& V0, Eigen::MatrixXi& T, double radius);


Rig* get_rig(std::string rig_file, Eigen::MatrixXd& V0, Eigen::MatrixXi& T, double radius, std::string weights, std::string weights_type);