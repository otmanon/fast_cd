#pragma once
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "igl/slice.h"
#include <Eigen/Dense>
#include "igl/volume.h"
#include "igl/massmatrix.h"
#include "interweaving_matrix.h"
#include "igl/repdiag.h"
#include <iostream>
/*
Computes the vector gradient operator K such that 
f = K(u + x), where f is the deformation gradient of the dispalcement map u

Inputs:
X - n x 3 vertex geometry
T - T x 4 tet mesh indicies

Output L
K - 9T x 3n tetrahedron gradient operator, following the matlab column ordering conveneiton (xxyyzz)
*/
void vector_gradient_operator(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,  Eigen::SparseMatrix<double>& K)
{

	int k = T.rows();
	int n = X.rows();



	
	K.resize(k * 9, n * 3);
	std::vector<Eigen::Triplet<double>> tripletList;

	//Eigen::Vector4i v_inds;
	Eigen::MatrixXd K_t, X_t, B_t;	//local for each tet

	B_t.resize(4, 3);
	B_t << -1, -1, -1,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	for (int t = 0; t < k; t++)
	{
		const Eigen::VectorXi v_inds = T.row(t); //indices of vertices belonging to this tet.
		igl::slice(X, v_inds, 1, X_t);
		X_t = X_t.transpose().eval();

		K_t = B_t * (X_t * B_t).inverse();

		for (int c = 0; c < K_t.cols(); c++)
		{
			Eigen::VectorXd K_t_c = K_t.col(0);
			for (int i = 0; i < T.cols(); i++)
			{
				int v_ind = T(t, i);

				tripletList.push_back(Eigen::Triplet<double>(3 * t + c * 3 * k + 0, v_ind, K_t(i, c)));
				tripletList.push_back(Eigen::Triplet<double>(3 * t + c * 3 * k + 1, v_ind + n, K_t(i, c)));
				tripletList.push_back(Eigen::Triplet<double>(3 * t + c * 3 * k + 2, v_ind + 2 * n, K_t(i, c)));

			}
			//each entry in K(r, c) contributes to 3 other places ... one for x, one for y, and one for z for the corresponding vertex i.
		}
	
	}
	K.setFromTriplets(tripletList.begin(), tripletList.end());

}