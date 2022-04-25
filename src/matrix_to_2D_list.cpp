#include "matrix_to_2D_list.h"

std::vector<std::vector<double>> matrix_to_2D_list(Eigen::MatrixXd& A)
{
	//list of rows
	std::vector<std::vector<double>> a;
	a.resize(A.rows());
	
	std::vector<double> row;
	row.resize(A.cols());
	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
		{
			row[j] = A(i, j);
		}
		a[i] = row;
	}
	return a;
}

