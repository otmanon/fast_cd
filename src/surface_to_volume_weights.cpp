#include "surface_to_volume_weights.h"
#include "prolongation.h"
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/procrustes.h>
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

	W = W.cwiseMax(0.0);
	W = W.cwiseMin(1.0);
	//normalize weight matrix
	Eigen::VectorXd sums = W.rowwise().sum();
	W.array().colwise() /= W.rowwise().sum().array();
	return W;
}


Eigen::MatrixXd surface_to_volume_weights(Eigen::MatrixXd& surfaceW, Eigen::MatrixXd& surfaceV, Eigen::MatrixXd& X, Eigen::MatrixXi& T)
{
	Eigen::MatrixXd W;
	Eigen::SparseMatrix<double> C, Aeq;
	igl::cotmatrix(X, T, C);
	C *= -1;
	Eigen::MatrixXd Beq;

	Eigen::MatrixXi F;
	Eigen::MatrixXd X_s, X_d, X_t;
	igl::boundary_facets(T, F);
	Eigen::VectorXi uniqueI;
	igl::unique(F, uniqueI);
	X_s = igl::slice(X, uniqueI, 1);
	Eigen::Matrix3d R;
	Eigen::VectorXd t;
	double scale;
	igl::procrustes(X_s, surfaceV, true, true, scale, R, t);
	X_d = ((scale * R * X_s.transpose()).colwise() + t).transpose();
	Eigen::VectorXi bI, all; 
	Eigen::VectorXd sqrD;
	Eigen::MatrixXd CP;

	X_t = ((scale * R * X.transpose()).colwise() + t).transpose();
	all = igl::colon<int>(0, X.rows() -1);
	//map surface vertices to closest volume vertices. These are now boundary conditions
	igl::point_mesh_squared_distance(surfaceV, X_t , all, sqrD, bI, CP);


	Eigen::MatrixXi BI = bI.replicate(1, surfaceW.rows());
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(C.rows(), surfaceW.cols());
	igl::min_quad_with_fixed_data<double> data;
	igl::min_quad_with_fixed_precompute(C, bI, Aeq, true, data);
	igl::min_quad_with_fixed_solve(data, B, surfaceW, Beq, W);

	W = W.cwiseMax(0.0);
	W = W.cwiseMin(1.0);
	//normalize weight matrix
	Eigen::VectorXd sums = W.rowwise().sum();
	W.array().colwise() /= W.rowwise().sum().array();
	return W;
}