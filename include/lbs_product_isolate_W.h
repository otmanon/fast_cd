#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <igl/blkdiag.h>
#include <igl/cat.h>
#include <igl/slice.h>
#include <igl/diag.h>
using namespace Eigen;
using namespace igl;
/*
	Consider a matrix multiplication like so C * B, where B is a linear blend skinning matrix. More precisle
	rewrite the LBS matrix C kron( I_3, (kron(1_w^T, [V 1]) .* kron(W, 1_4^T))) . 
	
	It is possible to rewrite this as a matrix multiplication acting directly on W
	
	CB = reshape(AW, j,  w*12).

	CB -> j x w * d * (d+1)
	This method returns the matrix A above. Follows derivation in overleaf
	and matlab 

	Input:
	V - > n x d matrix of rest pose vertex positions
	C - > j x m linear matrix applied on B
	Output
	A - > j*d*(d+ 1) x w   equivalent matrix acting on W
	*/
void lbs_product_isolate_W(const MatrixXd& V, const MatrixXd& C, MatrixXd& A)
{
	double d = V.cols();
	int n = V.rows();
	int m = C.cols();
	int r = C.rows();
	MatrixXd U = MatrixXd::Zero(n, d + 1);
	U.block(0, 0, n, d) = V;
	U.col(d) = VectorXd::Ones(n);
	
	MatrixXd G;

	for (int i = 0; i < d; i++)
	{
		MatrixXd Cd;
		
		for (int j = 0; j < d + 1; j++)
		{
			MatrixXd s = C.block(0, n * i, r,  n);
			blkdiag({ Cd, s }, Cd);
		}
		G = igl::cat(1, G, Cd);
	}

	SparseMatrix<double> U_tilde = SparseMatrix<double>(0, n);
	for (int i = 0; i < d + 1; i++)
	{
		SparseMatrix<double> Us;
		diag(U.col(i), Us);
		U_tilde = cat(1, U_tilde, Us);
	}
	//finally after all that work.
	A = G * U_tilde;
}