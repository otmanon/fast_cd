#pragma once
#include "null_rig.h"

NullRig::NullRig(Eigen::MatrixXd& X,
	Eigen::MatrixXi& T) : X(X), T(T)
{
	J.resize(X.rows() * 3, 0);
	J.setZero();
	p0.resize(0);
	N.resize(X.rows() * 3, X.rows() * 3);
	N.setIdentity();
};


void NullRig::get_rig_motion(Eigen::VectorXd& p, Eigen::VectorXd& ur)
{
	ur.resize(X.rows()*X.cols());
	ur.setZero();
}


void NullRig::get_rig_params(Eigen::VectorXd& p)
{
	p.resizeLike(p0);
}

