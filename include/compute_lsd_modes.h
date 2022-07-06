#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
/// <summary>
/// 
/// Computes modes using a similar linear subspace design methodology... scatter points around the mesh using farthest point sampling. number of points scattered is the number of modes.
/// Eeach mode corresponds to motion moving one point fully (w=1) and all other constrained points are zero (w=0).  Interpolate between thes econstraints by minimizing a laplacian^2 energy (aka biharmonic)
/// </summary>
/// <param name="V0"></param>
/// <param name="T"></param>
/// <param name="num_modes"></param>
/// <param name="B"></param>
void compute_lsd_modes(const Eigen::MatrixXd& V0, const  Eigen::MatrixXi& T, const Eigen::SparseMatrix<double>& J, int num_modes, Eigen::MatrixXd& B);