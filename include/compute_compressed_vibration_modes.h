#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
/*
Computes sparse vibration modes according to https://graphics.tudelft.nl/Publications-new/2017/BH17/comVibModes.pdf

Input:
H - hessian matrix
M - mass matrix
(Optional)
mu - paramter that ways sparsity-inducing penalty term
to_convergence - bool. If true, runs the algorithm until each mode converges.(default false)
max_iters - if to_convergence is false, the maximum number of iterations to spend on each mode
*/
void compute_compressed_vibration_modes(Eigen::SparseMatrix<double>& H, Eigen::SparseMatrix<double>& M, int num_modes, Eigen::MatrixXd& U,  double mu = 1e-2, bool to_convergence = false, int max_iters = 10);