#include "quadratic_rotated_matrix.h"
#include "constant_term_rotation_matrix.h"
QuadraticRotatedMatrix::QuadraticRotatedMatrix(const Eigen::MatrixXd& U, const Eigen::MatrixXd& V, const Eigen::SparseMatrix<double>& H, const Eigen::MatrixXd& W)
{
	int num_b = W.cols();
	UEHEV.resize(num_b);

	Eigen::SparseMatrix<double> E_uv, E_xy;
	for (int b = 0; b < num_b; b++)
	{
		UEHEV[b].resize(3);
		for (int u = 0; u < 3; u++)
		{
			UEHEV[b][u].resize(3);
			for (int v = 0; v < 3; v++)
			{
				UEHEV[b][u][v].resize(3);
				constant_term_rotation_matrix(W, u, v, b, E_uv);
				Eigen::MatrixXd D_temp = U * E_uv * H; 
				for (int x = 0; x < 3; x++)
				{
					UEHEV[b][u][v][x].resize(3);
					for (int y = 0; y < 3; y++)
					{
						constant_term_rotation_matrix(W, x, y, b, E_xy);
						const Eigen::MatrixXd A_temp = D_temp * E_xy * V; 
						UEHEV[b][u][v][x][y] = A_temp;
					}
				}
			}
		}
	}
}


void QuadraticRotatedMatrix::update(const Eigen::VectorXd& p, Eigen::MatrixXd& D)
{
	int num_b = p.rows() / 12;

	D.resize(UEHEV[0][0][0][0][0].rows(), UEHEV[0][0][0][0][0].cols());
	D.setZero();
	int xyi, uvi;
	for (int b = 0; b < num_b; b++)
	{
		for (int u = 0; u < 3; u++)
		{
			for (int v = 0; v < 3; v++)
			{
				uvi = 12 * b + 4 * u + v;
				for (int x = 0; x < 3; x++)
				{
					for (int y = 0; y < 3; y++)
					{
						xyi = 12 * b + 4 * x + y;
						D += (UEHEV[b][u][v][x][y] * p(xyi)* p(uvi));
					}
				}
			}
		}
	}
}

