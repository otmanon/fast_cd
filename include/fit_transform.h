#pragma once
#include "Eigen/Core"

/*
Given two versions of the same mesh, finds the scale mapping X to Y such that
SX = Y
*/
Eigen::Matrix3d fit_scale(Eigen::MatrixXd& X, Eigen::MatrixXd& Y)
{
	Eigen::RowVector3d s = (Y.colwise().maxCoeff() - Y.colwise().minCoeff()).array() / (X.colwise().maxCoeff() - X.colwise().minCoeff()).array();
	Eigen::Matrix3d S;
	S << s(0), 0, 0,
		0, s(1), 0,
		0, 0, s(2);

	return S;
}