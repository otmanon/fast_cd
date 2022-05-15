#include "get_rig_pinned_motion.h"
#include "matrix4f_from_parameters.h"
#include <igl/slice.h>
#include <igl/cat.h>
void get_rig_pinned_motion(Eigen::VectorXd& p_rest, Eigen::VectorXd& p_rel, Eigen::MatrixXd& X, std::vector<Eigen::VectorXi> bI, Eigen::VectorXd& bc)
{
	// used in satisfying S u = bc. bc must be ordered the same as the columns of S
	int ci = 0;
	Eigen::VectorXi I;
	Eigen::MatrixXd BC_1_tmp, BC_tmp, BC;
	Eigen::MatrixXd X_b, X_b1;
	Eigen::Matrix4d T;

	Eigen::Matrix4f A_inv;
	Eigen::Matrix4f A;

	for (int b = 0; b < bI.size(); b++)
	{
		//get transofrmation matrix
		Eigen::Matrix4f A_rest = matrix4f_from_parameters(p_rest, b);
		Eigen::Matrix4f A_inv = A_rest.inverse();
		Eigen::Matrix4f A = matrix4f_from_parameters(p_rel, b);
		Eigen::Matrix4d T = (A).cast<double>();

		I = bI[b];
		igl::slice(X, I, 1, X_b);

		Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(X_b.rows(), 1);
		igl::cat(2, X_b, ones, X_b1);

		BC_1_tmp = T * X_b1.transpose();
		BC_1_tmp.transposeInPlace();
		BC_tmp = BC_1_tmp.block(0, 0, BC_1_tmp.rows(), 3);

		BC.conservativeResize(BC.rows() + BC_tmp.rows(), 3);
		BC.bottomRows(BC_tmp.rows()) = BC_tmp - X_b;
	}

	bc = Eigen::Map<Eigen::VectorXd>(BC.data(), BC.rows() * 3);
}

