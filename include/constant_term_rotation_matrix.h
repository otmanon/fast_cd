#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
void constant_term_rotation_matrix(const Eigen::MatrixXd& W, const int u,const  int v, const int b, Eigen::SparseMatrix<double>& E_uv);