#pragma once
#include <Eigen\Core>

/*
Rig controller class for a null rig that simply sets and identifies motion of vertices centered at the C.O.M of X inside a ball of radius r
*/
class NullRigControllerScriptedConstrainedBall
{
public:
	NullRigControllerScriptedConstrainedBall(Eigen::MatrixXd& X, double r);

	void scripted_motion(int step, Eigen::MatrixXd& U);

	
public:
	double r;
	Eigen::VectorXi bI;	//constrained Indices
	Eigen::MatrixXd bc0;
};