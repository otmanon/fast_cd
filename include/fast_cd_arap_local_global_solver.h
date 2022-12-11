#pragma once
#include "cd_arap_local_global_solver.h"
#include "fast_cd_sim_params.h"
#include "fast_cd_arap_precomp.h"

struct fast_cd_arap_local_global_solver : cd_arap_local_global_solver
{

	LLT<MatrixXd> llt_solver;
	LDLT<MatrixXd> ldlt_solver;

	fast_cd_arap_local_global_solver() {};

	int prev_solve_iters;
	double prev_res;

	/*
	Constructs arap local global solver object used to solve 
	dynamics quickly for fast Complementary DYnamoics

	Inputs : 
		A - m x m system matrix
		Aeq - c x m constraint rows that enforece Aeq z = b 
				as linear equality constraints
		 p - cd_arap_local_global_solver_params
	*/
	fast_cd_arap_local_global_solver(const MatrixXd& A, const MatrixXd& Aeq, 
		const cd_arap_local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		MatrixXd H;
		augment_with_linear_constraints(A, Aeq, H);
		if (Aeq.rows() > 0)
			ldlt_solver.compute(H);
		else
			llt_solver.compute(H);

		prev_solve_iters = 0;
		prev_res = std::numeric_limits<double>::max();
	}


	/*
	Constructs arap local global solver object used to solve
	dynamics quickly for fast Complementary DYnamoics

	Inputs :
		A - m x m system matrix
		Aeq - c x m constraint rows that enforece Aeq z = b
				as linear equality constraints
		run_solver_to_convergence - (bool)
		max_iters - (int)
		convergence_threshold - (double)
						where to stop if ||res||2 drops below this
	*/
	fast_cd_arap_local_global_solver(const MatrixXd& A, const MatrixXd& Aeq,
		bool run_solver_to_convergence, int max_iters,
		double convergence_threshold)
	{
		this->p = cd_arap_local_global_solver_params(
			run_solver_to_convergence,
			max_iters, convergence_threshold);
		VectorXi bI;
		MatrixXd H;
		augment_with_linear_constraints(A, Aeq, H);
		if (Aeq.rows() > 0)
			ldlt_solver.compute(H);
		else
			llt_solver.compute(H);

		prev_solve_iters = 0;
		prev_res = std::numeric_limits<double>::max();
	}

	VectorXd solve(const VectorXd& z, const fast_cd_sim_params& params, const fast_cd_arap_dynamic_precomp& dp, const fast_cd_arap_static_precomp& sp)
	{
		Eigen::VectorXd z_next = z, z_prev = z;
		if (p.to_convergence)
		{
			int iter = 0;
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r = local_step(z_next, dp, sp);
				z_next = global_step(z_next, params, dp, sp, r);
				res = (z_next - z_prev).norm();
				iter += 1;
			} while (res > p.threshold);
			prev_res = res;
			prev_solve_iters = iter;
		}
		else
		{
			int iter = 0;
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r = local_step(z_next, dp, sp);
				z_next = global_step(z_next, params, dp, sp, r);
				res = (z_next - z_prev).norm();
				iter += 1;
				if (iter >= p.max_iters)
					break;
			} while (res > p.threshold);
			prev_res = res;
			prev_solve_iters = iter;
		}
		return z_next;
	};

	VectorXd local_step(const VectorXd& z, const fast_cd_arap_dynamic_precomp& dp, const fast_cd_arap_static_precomp& sp)
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

	VectorXd global_step(const VectorXd& z, const fast_cd_sim_params& params, const fast_cd_arap_dynamic_precomp& dp, const fast_cd_arap_static_precomp& sp, const VectorXd& r)
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