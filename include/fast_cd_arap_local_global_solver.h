#pragma once
#include "cd_arap_local_global_solver.h"
#include "fast_cd_sim_params.h"
#include "fast_cd_arap_precomp.h"

struct fast_cd_arap_local_global_solver
{
	cd_arap_local_global_solver_params p;
	VectorXd z;
	LLT<MatrixXd> llt_solver;
	LDLT<MatrixXd> ldlt_solver;

	fast_cd_arap_local_global_solver() {};
	fast_cd_arap_local_global_solver(MatrixXd& A, MatrixXd& Aeq, cd_arap_local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		MatrixXd H;
		augment_with_linear_constraints(A, Aeq, H);
		if (Aeq.rows() > 0)
			ldlt_solver.compute(H);
		else
			llt_solver.compute(H);
	}

	VectorXd solve(VectorXd& z, fast_cd_sim_params& params, fast_cd_arap_dynamic_precomp& dp, fast_cd_arap_static_precomp& sp)
	{
		Eigen::VectorXd z_next;
		for (int i = 0; i < p.max_iters; i++)
		{
			VectorXd r = local_step(z, dp, sp);
			z = global_step(z, params, dp, sp, r);
		}
		return z;
	};

	VectorXd local_step(VectorXd& z, fast_cd_arap_dynamic_precomp& dp, fast_cd_arap_static_precomp& sp)
	{
		VectorXd f = sp.GmKB * z + dp.GmKur + sp.GmKx;
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

	VectorXd global_step(VectorXd& z, fast_cd_sim_params& params, fast_cd_arap_dynamic_precomp& dp, fast_cd_arap_static_precomp& sp, VectorXd& r)
	{
		VectorXd inertia_grad = -dp.BMy;

		VectorXd arap_grad = dp.BCur + sp.BCx - sp.G1VKB.transpose() * r;

		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad + dp.f_ext;
		VectorXd Y;

		Eigen::VectorXd rhs = igl::cat(1, (- g).eval(), dp.bc);
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