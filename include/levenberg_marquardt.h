#pragma once
#include <Eigen/Core>
#include <functional>

// Conduct `max_iters` iterations of  "projected Levenberg-Marquardt descent
//
// Inputs:
//   f	     function that measures energy
//   df/da    function that computes gradient of our energy
//   dx/da		function that computes the kinematics jacobian
//   proj_z  function that projects z onto the set of feasible values
//   max_iters	
//	 tol
//   mu			--blending parameter that interpolated between gradient descent and gauss newton  
//   z  #z vector of initial z values
// Outputs
//   z  #z vector of optimized z values
void levenberg_marquardt(
	const std::function<double(const Eigen::VectorXd&)>& f,
	const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& dfda,
	const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dxda,			//
	const std::function<void(Eigen::VectorXd&)>& proj_z,
	const int max_iters,
	const double tol,
	const double mu,					//parameter that blends between gradient descent and gauss-newton
	Eigen::VectorXd& z);

