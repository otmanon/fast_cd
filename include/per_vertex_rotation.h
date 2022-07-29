#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/blkdiag.h>
#include "interweaving_matrix.h"
/*
Given rig parameters and rig weights, get a per vertex transformation matrix (may include shears), without the translational part.
*/
void per_vertex_rotation(const Eigen::VectorXd& p, const Eigen::MatrixXd& W, Eigen::SparseMatrix<double>& R, Eigen::MatrixXd& R_stack);