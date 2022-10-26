#pragma once


#include <Eigen/Core>
#include <Eigen/Sparse>

/// <summary>
/// Builds linear blend skinning jacobian matrix. J = kron( I3, [V 1] kron 1_w'  cdot kron kron 1_4' cdot
/// </summary>
/// <param name="V"></param>
/// <param name="W"></param>

using namespace Eigen;
void lbs_jacobian(const MatrixXd& V, const MatrixXd& W, SparseMatrix<double>& J);