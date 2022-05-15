#include "fast_complementary_dynamics_constrained_gevp.h"
#include "arap_hessian.h"
#include "momentum_leaking_matrix.h"
#include "sparse_diag.h"
#include "augment_with_linear_constraints.h"

#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/massmatrix.h>
#include "complementary_equality_constraint.h"
void fast_complementary_dynamics_constrained_geivp(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,
	const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& H_exp, Eigen::SparseMatrix<double>& M_exp)
{
    Eigen::SparseMatrix<double> H, M;
	//arap hessian
	arap_hessian(X, T, H);
	H *= -1;


    //augment hessian with complementary constraint jacobians.
    Eigen::SparseMatrix<double> A_eq;

    complementary_equality_constraint(X, T, J, A_eq);
    augment_with_linear_constraints(H, A_eq, H_exp);

	
    //append tiny diagonal entries to sparse matrix...
    igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    M = igl::repdiag(M, 3);
    Eigen::VectorXd m = M.diagonal();
    Eigen::VectorXd z = Eigen::VectorXd::Zero(J.cols());
    Eigen::VectorXd diag = igl::cat(1, m, z);
    M_exp = sparse_diag(diag);


    Eigen::SparseMatrix<double> I;
    I.resize(H.rows(), H.rows());
    I.setIdentity();

   // H += 1e-8 * I;

}
