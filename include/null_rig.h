#pragma once
#include "rig.h"
class NullRig : public Rig
{
public:
	NullRig(Eigen::MatrixXd& X,
	Eigen::MatrixXi& T);

	
	//TODO: distinguish between get_rig_params. Get rig motion goes from the input params and should output a full dimensional motion
	void get_rig_motion(Eigen::VectorXd& p, Eigen::VectorXd& ur);
	void get_rig_params(Eigen::VectorXd& p);

public:

	//Rest pose geometry
	Eigen::MatrixXd X;
	Eigen::MatrixXi T;


	//Jacobian
	Eigen::SparseMatrix<double> J;

	//rig parameters, and original rig parameters. (might not even have to keep track of this)
	Eigen::VectorXd p0;
};