#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "igl/slice.h"
#include <Eigen/Dense>
#include "igl/volume.h"
#include "igl/massmatrix.h"
#include "interweaving_matrix.h"
#include "igl/repdiag.h"
#include <iostream>
void deformation_gradient_from_u_prefactorized_matrices(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, Eigen::VectorXd& H, Eigen::SparseMatrix<double>& K, Eigen::SparseMatrix<double>& M)
{

	int k = T.rows();
	int n = X.rows();

	Eigen::VectorXd vol;
	igl::volume(X, T, vol);

	H.resize(k * 9);
	K.resize(k * 9, n * 3);
	std::vector<Eigen::Triplet<double>> tripletList;

	//Eigen::Vector4i v_inds;
	Eigen::MatrixXd K_t, X_t, B_t;	//local for each tet
	Eigen::Matrix3d H_t = Eigen::Matrix3d::Identity();	//local for each tet
	double vol_t;
	B_t.resize(4, 3);
	B_t << -1, -1, -1,
			1,  0,  0,
			0,  1,  0,
			0,  0,  1;


	for (int t = 0; t < k; t++)
	{
		vol_t = vol(t);
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

				tripletList.push_back(Eigen::Triplet<double>(3 * t + c*3*k + 0, v_ind      , K_t(i, c)));
				tripletList.push_back(Eigen::Triplet<double>(3 * t + c*3*k + 1, v_ind +   n, K_t(i, c)));
				tripletList.push_back(Eigen::Triplet<double>(3 * t + c*3*k + 2, v_ind + 2*n, K_t(i, c)));

			}
					//each entry in K(r, c) contributes to 3 other places ... one for x, one for y, and one for z for the corresponding vertex i.
		}
		//Push entries in H
		H.middleRows(3 * t          , 3) = H_t.col(0);
		H.middleRows(3 * t + 3*k    , 3) = H_t.col(1);
		H.middleRows(3 * t + 3*k * 2, 3) = H_t.col(2);


	}
	K.setFromTriplets(tripletList.begin(), tripletList.end());

	//igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	M.resize(T.rows(), T.rows());
	M.setIdentity();
	M.diagonal() = vol;
	Eigen::SparseMatrix<double> S; // interweaving matrix

	interweaving_matrix(M.rows(), 3, S);

	M = igl::repdiag(M, 3);

	M =  S.transpose() * M * S ;

	M = igl::repdiag(M, 3);


	std::vector<Eigen::Triplet<double>> mlist;
	Eigen::SparseMatrix<double> M2;
	
	for (int t = 0; t < T.rows(); t++)
	{
		for (int i = 0; i < 3; i++)
		{
		//	double entry1, entry2, entry3, entry4;
		//	entry1 = M.coeffRef(3 * t + i, 3 * t + i);
		//	entry4 = M.coeffRef(3 * t + i + 1, 3 * t + i + 1);
		//	entry2 = M.coeffRef(3 * t + T.rows() + i, 3 * t + T.rows() + i);
		//	entry3 = M.coeffRef(3 * t + 2*T.rows() + i, 3 * t + 2*T.rows() + i);
			mlist.push_back(Eigen::Triplet<double>(3 * t +              i, 3 * t +              i, vol(t)));
			mlist.push_back(Eigen::Triplet<double>(3 * t +   3*T.rows() + i, 3 * t + 3*T.rows() +   i, vol(t)));
			mlist.push_back(Eigen::Triplet<double>(3 * t + 6*T.rows() + i, 3 * t + 6*T.rows() + i, vol(t)));
		}
	}
	M2.resize(T.rows()*9, T.rows()*9);
	M2.setFromTriplets(mlist.begin(), mlist.end());
	//std::cout << "Mass matrix error : " << (M2 - M).squaredNorm() << std::endl;
	M = M2;
}