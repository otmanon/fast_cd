#pragma once

#include "surface_to_volume_weights.h"
#include "rig_parameters.h"
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/procrustes.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
#include <Eigen/Core>
using namespace Eigen;



/*
Given a surface rig with parameters P0, fit it to the mesh defined by V, T. 
The new rig bones might be changing scale, rotation and centroid....
Inputs :
V - mesh geometry
T - tet indices
Vs - intiial surface geometry associated with rig
P0 - initial rig parameters
lengths - intiial rig lengths
Outpus:
P0 - fitted rig parameters
*/
void fit_rig_to_mesh(const MatrixXd& V, const MatrixXi& T, const MatrixXd& Vs, MatrixXd& P0)
{
	Eigen::MatrixXi F;

	Eigen::MatrixXd Xs;
	igl::boundary_facets(T, F);
	Eigen::VectorXi uniqueI;
	igl::unique(F, uniqueI);
	igl::slice(V, uniqueI, 1, Xs);

	double ss; Matrix3d R; Vector3d t;
	igl::procrustes(Xs, Vs, true, true, ss, R, t);


	Eigen::Matrix4d A = Matrix4d::Identity();
	A.block(0, 0, 3, 3) =  ss*R.transpose();
	A.block(0, 3, 3, 1) = t;
	A = A.inverse().eval();
	MatrixXd B = A.block(0, 0, 3, 4);
	transform_rig_parameters(P0, B);

	//lengths *= ss;
}

/*
Given a surface rig with parameters P0, fit it to the mesh defined by V, T.
The new rig bones might be changing scale, rotationand centroid....
Inputs :
V - mesh geometry
T - tet indices
Vs - intiial surface geometry associated with rig
P0 - initial rig parameters
lengths - intiial rig lengths
Outpus :
P0 - fitted rig parameters
B - the affine matrix that fits the rig to the mesh
*/
void fit_rig_to_mesh(const MatrixXd& V, const MatrixXi& T, const MatrixXd& Vs, MatrixXd& P0, MatrixXd& B)
{
	Eigen::MatrixXi F;

	Eigen::MatrixXd Xs;
	igl::boundary_facets(T, F);
	Eigen::VectorXi uniqueI;
	igl::unique(F, uniqueI);
	igl::slice(V, uniqueI, 1, Xs);

	double ss; Matrix3d R; Vector3d t;
	igl::procrustes(Xs, Vs, true, false, ss, R, t);


	Matrix4d A = Matrix4d::Identity();
	A.block(0, 0, 3, 3) = ss * R.transpose();
	A.block(0, 3, 3, 1) = t;
	A = A.inverse().eval();
	B = A.block(0, 0, 3, 4);
	transform_rig_parameters(P0, B);

}


void fit_rig_to_mesh_vertices(const MatrixXd& V, const MatrixXd& Vs, MatrixXd& P0, MatrixXd& B)
{
	Eigen::MatrixXi F;


	double ss; Matrix3d R; Vector3d t;
	igl::procrustes(V, Vs, true, false, ss, R, t);

	Matrix4d A = Matrix4d::Identity();
	A.block(0, 0, 3, 3) = ss * R.transpose();
	A.block(0, 3, 3, 1) = t;
	A = A.inverse().eval();
	 B = A.block(0, 0, 3, 4);
	transform_rig_parameters(P0, B);


}
