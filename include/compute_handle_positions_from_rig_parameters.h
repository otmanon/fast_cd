#pragma once
#include <Eigen\Core>

//these are row order flattened
void compute_handle_positions_from_parameters(Eigen::VectorXd& p, Eigen::MatrixXd& V);