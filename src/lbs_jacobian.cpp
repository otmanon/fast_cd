#include "lbs_jacobian.h"
#include "kron.h"


#include <igl/cat.h>
#include <igl/repdiag.h>

#include <assert.h>

/// 
/// Builds linear blend skinning jacobian matrix. J = kron( I3, [V 1] kron 1_w'  cdot kron kron 1_4' cdot
/// Inputs
/// V -> n x 3 geometry, positions at rest state
/// W -> n x b bones weights
/// J -> 3n x 12b linear blend skinning jacobian
///
using namespace Eigen;

void lbs_jacobian(const MatrixXd& V, const MatrixXd& W, MatrixXd& J)
{

	assert(V.rows() == W.rows() && "Weights should have same number of rows as W!");

	assert(V.cols() == 3 && "Can only handle 3D case for now");

	int n = V.rows();
	int dim = V.cols();
	int b = W.cols();
	MatrixXd ones_n = MatrixXd::Ones(n, 1);
	MatrixXd U = igl::cat(2, V, ones_n);

	MatrixXd U_exp;
	MatrixXd ones_b = MatrixXd::Ones(b, 1).transpose(); //row vector of b ones
	kron(ones_b, U, U_exp);

	MatrixXd W_exp;
	MatrixXd ones_dp1 = MatrixXd::Ones(dim + 1, 1).transpose(); //row vector of d+1 ones
	kron(W, ones_dp1, W_exp);

	Eigen::MatrixXd J_compact = U_exp.array() * W_exp.array();  //component wise product

	J = igl::repdiag(J_compact, dim);

}


void lbs_jacobian(const MatrixXd& V, const MatrixXd& W, SparseMatrix<double>& J)
{
	MatrixXd Jd;
	lbs_jacobian(V, W, Jd);
	J = Jd.sparseView(); // TODO: have a separate method for more quickly building the sparse matrix J
}