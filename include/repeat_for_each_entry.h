#pragma once
#include "kron.h"
#include <Eigen/Core>

using namespace Eigen;
/*
takes each entry in a and repeats it n, times, and tiles it m times
*/
VectorXd  repeat_for_each_entry(VectorXd& a, int n, int m)
{
	MatrixXd a_tmp1, a_tmp2;
	VectorXd a_exp;
	VectorXd ones1 = VectorXd::Ones(n);
	VectorXd ones2 = VectorXd::Ones(m);
	kron(a, ones1, a_tmp1);
	kron(ones2, a_tmp1, a_tmp2);
	a_exp = a_tmp2;
	return a_exp;
};

VectorXi  repeat_for_each_entry(VectorXi& a, int n, int m)
{
	VectorXd ad = a.cast<double>();
	VectorXd ad_exp = repeat_for_each_entry(ad, n, m);
	VectorXi a_exp = ad_exp.cast<int>();
	return a_exp;
};