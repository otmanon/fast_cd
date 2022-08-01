#include "farthest_point_sampling.h"
#include <igl/heat_geodesics.h>
#include <igl/slice.h>
void farthest_point_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const int num_samples, Eigen::VectorXi& I, Eigen::MatrixXd& P)
{
	igl::HeatGeodesicsData<double> data;
	igl::heat_geodesics_precompute(V, T, data);

	I.resize(1);
	
	Eigen::VectorXd D; //distances
	D = V.col(1);
	int ind;
	double max_d = D.maxCoeff(&ind);
	I(0) = ind; //initialize first point with highest point... why not!

	for (int i = 0; i < num_samples -1; i++)
	{
		igl::heat_geodesics_solve(data, I, D);

		double max_d = D.maxCoeff(&ind);
		I.conservativeResize(I.rows() + 1);
		I(I.rows() - 1) = ind;
	}

	//also get the vertex positions indexed by I
	igl::slice(V, I, 1, P);
}