#include "per_vertex_rotation.h"

void per_vertex_rotation(const Eigen::VectorXd& p, const Eigen::MatrixXd& W, Eigen::SparseMatrix<double>& R, Eigen::MatrixXd& R_stack)
{
	int num_b = p.rows() / 12;

	Eigen::SparseMatrix<double> S;
	interweaving_matrix(W.rows(), 3, S);

	R_stack.resize(3 * W.rows(), 3);
	R_stack.setZero();

	std::vector<Eigen::Triplet<double>> triplet_list;
	Eigen::Matrix3d A;
	for (int i = 0; i < W.rows(); i++)
	{
		for (int b = 0; b < num_b; b++)
		{
			const double w = W(i, b);
			A << p(12 * b + 0), p(12 * b + 1), p(12 * b + 2),
				p(12 * b + 4), p(12 * b + 5), p(12 * b + 6),
				p(12 * b + 8), p(12 * b + 9), p(12 * b + 10);
			A *= w;
			R_stack.block(3 * i, 0, 3, 3) += A;

			triplet_list.emplace_back(3 * i + 0, 3 * i + 0, w * p(12 * b + 0));
			triplet_list.emplace_back(3 * i + 0, 3 * i + 1, w * p(12 * b + 1));
			triplet_list.emplace_back(3 * i + 0, 3 * i + 2, w * p(12 * b + 2));


			triplet_list.emplace_back(3 * i + 1, 3 * i + 0, w * p(12 * b + 4));
			triplet_list.emplace_back(3 * i + 1, 3 * i + 1, w * p(12 * b + 5));
			triplet_list.emplace_back(3 * i + 1, 3 * i + 2, w * p(12 * b + 6));

			triplet_list.emplace_back(3 * i + 2, 3 * i + 0, w * p(12 * b + 8));
			triplet_list.emplace_back(3 * i + 2, 3 * i + 1, w * p(12 * b + 9));
			triplet_list.emplace_back(3 * i + 2, 3 * i + 2, w * p(12 * b + 10));
		}
	}

	R.resize(3 * W.rows(), 3 * W.rows());
	R.setFromTriplets(triplet_list.begin(), triplet_list.end());

	R = S * R * S.transpose();

	R_stack = S * R_stack;
}