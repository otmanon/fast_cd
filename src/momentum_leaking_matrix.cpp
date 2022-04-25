#include "momentum_leaking_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice_into.h>
#include <Eigen/Sparse>
#include "sparse_diag.h"
using namespace Eigen;
void momentum_leaking_matrix(const MatrixXd& V, const MatrixXi& F, fast_cd::MOMENTUM_LEAK_MATRIX method, SparseMatrix<double>& D )
{
	double dt = 1e-6;
	SparseMatrix<double> C, M;
	igl::cotmatrix(V, F, C);
	C *= -1;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	

	VectorXd Z0 = VectorXd::Zero(V.rows());
	VectorXi bI; MatrixXi bF; 
	igl::boundary_facets(F, bF);
	igl::unique(bF, bI);
	VectorXd Zb = VectorXd::Ones(bI.rows());
	igl::slice_into(Zb, bI, Z0);

	Eigen::SimplicialLLT<SparseMatrix<double>> solver;
	
	solver.compute(C + dt * M);
	Eigen::VectorXd Z = solver.solve(M * Z0);

	//recenter to lie between 0 and 1
	Z = Z.array() - Z.minCoeff();
	Z = Z.array() / Z.maxCoeff();

	//default is to have this be inverted
	if (method == fast_cd::MOMENTUM_LEAK_DIFFUSION)
	{
		Z = 1 - Z.array();
	}

	D = sparse_diag(Z);

}