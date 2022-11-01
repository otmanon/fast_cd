#pragma once

#include "surface_to_volume_weights.h"
#include "rig_parameters.h"
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/procrustes.h>
#include <igl/point_mesh_squared_distance.h>

#include <Eigen/Core>
using namespace Eigen;
/*
Fits rig parameters, surface weights to be well matched to the volume mesh X, T. 
This involves a few things:
	1. Mapping weights from a surface mesh to a volume mesh
	2. Scaling, rotating  and translating the rig to be well matched to the mesh X, T

Input :
W -> sn x b set of weights. 
P0 -> (dim+1)b x dim stacked rest pose rig parameters
V0 -> sn x 3 vertex positions matching the rig
X -> #V x dim tet mesh vertices
T -> #T x simplex_size tet mesh indices.

Output:
W -> Overwrwitten #V x b set of weights, across the entire tet mesh
P0 -> Overwritten (dim+1)b x dim rest pose rig parameters that have been scaled/ rotated/ fitted to match X, T
ss -> scalar Scale, done to fit geometry
R -> 3x3  rotation matrix used to fit geometry
t -> translation used to fit geometry

		If sn is == X.rows(), then nothing is done. If it isn't then the geometry of V0 is fitted to X, T and we solve for the  missing
							weights using a constrained diffusion optimization


*/
void fit_rig_to_mesh(MatrixXd& W, MatrixXd& P0, MatrixXd& V0, MatrixXd& X, MatrixXi& T, double& ss, Matrix3d& R, Vector3d& t)
{
	
	Eigen::MatrixXi F;
	
	Eigen::MatrixXd X_s, X_d, X_t;

	igl::boundary_facets(T, F);
	Eigen::VectorXi uniqueI;
	igl::unique(F, uniqueI);
	X_s = igl::slice(X, uniqueI, 1);

	
	igl::procrustes(X_s, V0, true, false, ss, R, t);
	X_d = ((ss* R * X_s.transpose()).colwise() + t).transpose();
	Eigen::VectorXi bI, all;
	Eigen::VectorXd sqrD;
	Eigen::MatrixXd CP;

	X_t = ((ss* R * X.transpose()).colwise() + t).transpose();
	all = igl::colon<int>(0, X.rows() - 1);
	//map surface vertices to closest volume vertices. These are now boundary conditions
	igl::point_mesh_squared_distance(V0, X_t, all, sqrD, bI, CP);

	W = surface_to_volume_weights(W, bI, X, T);

	Eigen::Matrix4d A = Matrix4d::Identity();
	A.block(0, 0, 3, 3) = ss * R.transpose();
	A.block(0, 3, 3, 1) = t;
	A = A.inverse().eval(); 
	MatrixXd B = A.block(0, 0, 3, 4);
	transform_rig_parameters(P0, B);
//	translate_rig_parameters(P0, t);
//	rotate_rig_parameters(P0, R);

}


void fit_rig_to_mesh(MatrixXd& W, MatrixXd& P0, MatrixXd& V0, MatrixXd& X, MatrixXi& T)
{
	Matrix3d R;
	double s;
	Vector3d t;
	fit_rig_to_mesh(W, P0, V0, X, T, s, R, t);

}

/*
Wrapper that exposes the affine transfomr matrix A
*/
void fit_rig_to_mesh(MatrixXd& W, MatrixXd& P0, MatrixXd& V0, MatrixXd& X, MatrixXi& T, MatrixXd& A)
{
	Matrix3d R;
	double s;
	Vector3d t;
	fit_rig_to_mesh(W, P0, V0, X, T, s, R, t);
	
	A = Matrix4d::Identity();
	A.block(0, 0, 3, 3) = s * R;
	A.block(0, 3, 3, 1) = t;
	A = A.inverse().eval();
	A = (A.block(0, 0, 3, 4)).eval();
}