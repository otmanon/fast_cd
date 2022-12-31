#pragma once
#include <Eigen/Sparse>


using namespace Eigen;
/*
Returns a matrix that when multiplied against a stack of flattened matrices,
returns the transposed stack of flattened matrices

n - (int) number of rows of each of the origin flattened matrix
m - (int) number of coloumns of each of the original unflattened matrix
d - (int) number of stacked matrices.
*/

void vectorized_transpose(int n, int m, int d, SparseMatrix<double>& AT)
{
	std::vector<Triplet<int>> tripletList;
	tripletList.reserve(n * m * d);
	int j = 0;
	
	//Assuming column order matlab-style flattening... outer loop are the columns
	for (int mi = 0; mi < m; mi++)
	{
		for (int di = 0; di < d; di++)
		{
			for (int ni = 0; ni < n; ni++)
			{

				int j2 = mi * (d * n) + di*n +  ni;
				int i = ni * (d * n) + di * n + mi; //inverse of j2
				tripletList.emplace_back(i, j, 1);
				j += 1;
			}
		}
	}
	AT.resize(n * m * d, n * m * d);
	AT.setFromTriplets(tripletList.begin(), tripletList.end());
	

	for (int trial = 0; trial < 100; trial++)
	{
		/// Test this out
		MatrixXd R, T;
		R.resize(d * n, m);
		T = R;
		for (int di = 0; di < d; di++)
		{
			Matrix3d v = Matrix3d::Random();
			R.block(di * n, 0, n, m) = v;
			T.block(di * n, 0, n, m) = v.transpose();
		}

		VectorXd r = Map<VectorXd>(R.data(), R.cols() * R.rows());
		VectorXd t = Map<VectorXd>(R.data(), T.cols() * T.rows());

		/*
		Make sure this is true:
		*/
		printf("Error : %d \n", (r - AT * t).norm());
	}
	

}