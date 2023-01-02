#pragma once
#include <Eigen/Core>
/*
From rest pose parameters p0, and deformed pose parameters sol_p, get p_rel
*/
void get_relative_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& p, Eigen::VectorXd& p_rel);