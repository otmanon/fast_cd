#include "list_2D_to_matrix.h"


Eigen::MatrixXd list_2D_to_matrix(std::vector<std::vector<double>> A_list)
{

	//list of rows
	Eigen::MatrixXd A;
	A.resize(A_list.size(), A_list[0].size());
	for (int i = 0; i < A.rows(); i++)
	{
		A.row(i) = Eigen::Map<Eigen::RowVectorXd>(&A_list[i][0], A.cols());
	}

	return A;
}