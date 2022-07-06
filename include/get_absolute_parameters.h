#pragma once
#include <Eigen/Core>
#include "matrix4f_from_parameters.h"
#include "update_parameters_at_handle.h"

/*
Given rest pose params and relative params, computes absolute params
*/
void get_absolute_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& p_rel, Eigen::VectorXd& p);