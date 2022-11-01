#pragma once
#include "cd_sim_params.h"
#include "cd_arap_precomp.h"
#include "augment_with_linear_constraints.h"

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <igl/cat.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
using namespace Eigen;
using namespace igl;

struct cd_arap_local_global_solver_params
{
	int max_iters;
	bool to_convergence;
	double threshold;

	cd_arap_local_global_solver_params() {};
	cd_arap_local_global_solver_params(bool to_convergence, int max_iters, double threshold)
	{
		this->to_convergence = to_convergence;
		this->max_iters = max_iters;
		this->threshold = threshold;
	}
};

struct cd_arap_local_global_solver
{
	cd_arap_local_global_solver_params p;
	VectorXd z;

	min_quad_with_fixed_data<double> data;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_solver;

	cd_arap_local_global_solver() {};
	cd_arap_local_global_solver(SparseMatrix<double>& A, SparseMatrix<double>& Aeq, cd_arap_local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		SparseMatrix<double> H;
		augment_with_linear_constraints(A, Aeq, H);
		ldlt_solver.compute(H);
		min_quad_with_fixed_precompute(A, bI, Aeq, false, data);
	}

	VectorXd solve(VectorXd& z, cd_sim_params& params, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		Eigen::VectorXd z_next;
		for (int i = 0; i < p.max_iters; i++)
		{
			VectorXd r = local_step(z, dp, sp);
			z = global_step(z, params, dp, sp, r);
		}
		return z;
	};

	VectorXd local_step(VectorXd& z, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		VectorXd f = sp.K * z + dp.Kur + sp.Kx;

		VectorXd f2 = sp.VK * z + dp.VKur + sp.VKx;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd F2_stack = Map<MatrixXd>(f2.data(), nt*3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		MatrixXd R2 = R;
		
		Matrix3d F, F2, rot, rot2;
		for (int i = 0; i < nt; i++)
		{
			F =  F_stack.block(3 * i, 0, 3, 3);// *sp.tet_vols(i);
		
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		return r;
	}

	VectorXd global_step(VectorXd& z, cd_sim_params& params, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp, VectorXd& r)
	{
		VectorXd inertia_grad = -dp.My;
		VectorXd arap_grad = dp.Cur + sp.Cx - sp.VK.transpose() * r;

		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad + dp.f_ext;
		VectorXd Y;

		Eigen::VectorXd rhs = -igl::cat(1, g, dp.bc);
		VectorXd z_next = ldlt_solver.solve(rhs);
		
		return z_next.topRows(params.X.rows()*3);
	}

};