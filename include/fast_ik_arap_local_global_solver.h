#pragma once
#include "fast_cd_arap_local_global_solver.h"

using namespace Eigen;
struct fast_ik_arap_local_global_solver 
	: fast_cd_arap_local_global_solver
{
	fast_ik_arap_local_global_solver(const MatrixXd& A, const MatrixXd& Aeq,
		const cd_arap_local_global_solver_params& p) :
		fast_cd_arap_local_global_solver(A, Aeq, p)
	{};


	VectorXd local_step(const VectorXd& z, const fast_cd_arap_dynamic_precomp& dp, const fast_cd_arap_static_precomp& sp)
	{
		VectorXd f = sp.GmKB * z;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		//	cout << F_stack << endl;
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());

		Matrix3d F, rot;
		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		return r;
	}

	VectorXd global_step(const VectorXd& z, const fast_cd_sim_params& params, const fast_cd_arap_dynamic_precomp& dp, const fast_cd_arap_static_precomp& sp, const VectorXd& r)
	{
		
		VectorXd x = Map<const VectorXd>(params.X.data(), params.X.rows() * params.X.cols());
		VectorXd arap_grad = -sp.G1VKB.transpose() * r;

		VectorXd g = arap_grad + dp.f_ext;
		VectorXd Y;

		Eigen::VectorXd rhs = igl::cat(1, (-g).eval(), dp.bc);
		VectorXd z_next;
		if (dp.bc.rows() > 0)
			z_next = ldlt_solver.solve(rhs);
		else
			z_next = llt_solver.solve(rhs);
		//min_quad_with_fixed_solve(data, g, VectorXd(), dp.bc, z_next);
		//min_quad_with_fixed<double>(sp.BAB, g, VectorXi(), VectorXd(), sp.AeqB, dp.bc);
		//arap = dpre.Lur + spre.Lx - spre.MK' * r;
		return z_next.topRows(params.B.cols());
	}

};