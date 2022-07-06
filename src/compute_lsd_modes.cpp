#include "compute_lsd_modes.h"
#include "complementary_equality_constraint.h"

#include <igl/invert_diag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
void compute_lsd_modes(const Eigen::MatrixXd& V0, const  Eigen::MatrixXi& T, const Eigen::SparseMatrix<double>& J, int num_modes, Eigen::MatrixXd& B)
{
	igl::min_quad_with_fixed_data<double> data;

	Eigen::SparseMatrix<double> M, C, Mi, A, Aeq, JTMD;

	igl::cotmatrix(V0, T, C);
	igl::massmatrix(V0, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	igl::invert_diag(M, Mi);

	A = C.transpose() * Mi * C;

	//need to fill Aeq... bottom will be with rig jacobianm  
	complementary_equality_constraint(V0, T, J, JTMD);

	Eigen::VectorXi _n;
	igl::min_quad_with_fixed_precompute(A, _n, Aeq, false, data);

}