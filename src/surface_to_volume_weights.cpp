#include "surface_to_volume_weights.h"
#include <Eigen/Sparse>
#include "igl/min_quad_with_fixed.h"
#include "igl/cotmatrix.h"
Eigen::MatrixXd surface_to_volume_weights(Eigen::MatrixXd& surfaceW, Eigen::VectorXi& bI, Eigen::MatrixXd& X, Eigen::MatrixXi& T)
{
	Eigen::MatrixXd W;
	Eigen::SparseMatrix<double> C, Aeq;
	igl::cotmatrix(X, T, C);
	C *= -1;
	Eigen::MatrixXd Beq;

 	Eigen::MatrixXi BI = bI.replicate(1, surfaceW.rows());
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(C.rows(), surfaceW.cols());
	igl::min_quad_with_fixed_data<double> data;
	igl::min_quad_with_fixed_precompute(C, bI, Aeq, true, data);
	igl::min_quad_with_fixed_solve(data, B, surfaceW, Beq, W);
	//normalize weight matrix
	Eigen::VectorXd sums = W.rowwise().sum();
	W.array().colwise() /= W.rowwise().sum().array();
	return W;
}