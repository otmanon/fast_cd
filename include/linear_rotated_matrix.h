#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

class LinearRotatedMatrix
{
public:
	//class that aids in the computation of D = U * R * V, where R is a 3nx3n block diagonal matrix of per vertex transformations. 
	//Key in this method is that we assume the per vertex transformations are linear functions of the rig parameters, as weighed by the rig Weights. 

	//This method constructs n prefactorized matrices that we then sum at run-time to quickly evaluate the per vertex quantity as O(rig parameter) cost. This method is a modification of the Fast Sandwich transform
	//where we don't assume multiple vertices have the same rotation matrix, but instead the rotation matrices are all different for each vertex, they are just a linear combination of rig parameters.
	/// <summary>
	/// Constructor for the fast sandiwch transformer
	/// </summary>
	/// <param name="U"> mx3n left hand matrix</param>
	/// <param name="V"> 3nxl right hand matrix</param>
	/// <param name="W"> rig weights</param>
	LinearRotatedMatrix(const Eigen::MatrixXd& U, const Eigen::MatrixXd& V, const Eigen::MatrixXd& W);
	LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::SparseMatrix<double>& U, const Eigen::MatrixXd& V, const Eigen::MatrixXd& W);

	LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::MatrixXd& U, const Eigen::SparseMatrix<double>& V, const Eigen::MatrixXd& W);

	LinearRotatedMatrix::LinearRotatedMatrix(const Eigen::SparseMatrix<double>& U, const Eigen::SparseMatrix<double>& V, const Eigen::MatrixXd& W);
	LinearRotatedMatrix() {};

	void update(const Eigen::VectorXd& p, Eigen::MatrixXd& D);

		
public:
	

	std::vector<std::vector<std::vector<Eigen::MatrixXd>>> UEV; // list of 9 matrices per each bone such that UEV_uvb p_uvb = U * T * V
};