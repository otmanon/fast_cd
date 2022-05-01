#pragma once
#include <Eigen/Core>
/*
Given bone parameters, and a skeleton hierarchy pI, get bone lengths
*/
double get_bone_length_scales(Eigen::VectorXd& p, Eigen::VectorXi& pI);