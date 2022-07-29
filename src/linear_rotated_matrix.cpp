#include "linear_rotated_matrix.h"
#include "interweaving_matrix.h"
#include "constant_term_rotation_matrix.h"


LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::MatrixXd& U, const Eigen::MatrixXd& V, const Eigen::MatrixXd& W)
{
	int num_b = W.cols();
	UEV.resize(num_b);

	Eigen::SparseMatrix<double> E_uv, E_xy;
	for (int b = 0; b < num_b; b++)
	{
		UEV[b].resize(3);
		for (int u = 0; u < 3; u++)
		{
			UEV[b][u].resize(3);
			for (int v = 0; v < 3; v++)
			{
				constant_term_rotation_matrix(W, u, v, b, E_uv);
			    UEV[b][u][v] = U * E_uv * V;
				
			}
		}
	}
}

LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::SparseMatrix<double>& U, const Eigen::MatrixXd& V, const Eigen::MatrixXd& W)
{
	int num_b = W.cols();
	UEV.resize(num_b);

	Eigen::SparseMatrix<double> E_uv, E_xy;
	for (int b = 0; b < num_b; b++)
	{
		UEV[b].resize(3);
		for (int u = 0; u < 3; u++)
		{
			UEV[b][u].resize(3);
			for (int v = 0; v < 3; v++)
			{
				constant_term_rotation_matrix(W, u, v, b, E_uv);
				UEV[b][u][v] = U * E_uv * V;

			}
		}
	}
}

LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::MatrixXd& U, const Eigen::SparseMatrix<double>& V, const Eigen::MatrixXd& W)
{
	int num_b = W.cols();
	UEV.resize(num_b);

	Eigen::SparseMatrix<double> E_uv, E_xy;
	for (int b = 0; b < num_b; b++)
	{
		UEV[b].resize(3);
		for (int u = 0; u < 3; u++)
		{
			UEV[b][u].resize(3);
			for (int v = 0; v < 3; v++)
			{
				constant_term_rotation_matrix(W, u, v, b, E_uv);
				UEV[b][u][v] = U * E_uv * V;

			}
		}
	}
}

LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::SparseMatrix<double>& U, const Eigen::SparseMatrix<double>& V, const Eigen::MatrixXd& W)
{
	int num_b = W.cols();
	UEV.resize(num_b);

	Eigen::SparseMatrix<double> E_uv, E_xy;
	for (int b = 0; b < num_b; b++)
	{
		UEV[b].resize(3);
		for (int u = 0; u < 3; u++)
		{
			UEV[b][u].resize(3);
			for (int v = 0; v < 3; v++)
			{
				constant_term_rotation_matrix(W, u, v, b, E_uv);
				UEV[b][u][v] = U * E_uv * V;

			}
		}
	}
}

void LinearRotatedMatrix::update(const Eigen::VectorXd& p, Eigen::MatrixXd& D)
{
	int num_b = p.rows() / 12;

	D.resize(UEV[0][0][0].rows(), UEV[0][0][0].cols());
	D.setZero();
	int xyi, uvi;
	for (int b = 0; b < num_b; b++)
	{
		for (int u = 0; u < 3; u++)
		{
			for (int v = 0; v < 3; v++)
			{
				uvi = 12 * b + 4 * u + v;
				D += UEV[b][u][v] * p(uvi);
			}
		}
	}
}

