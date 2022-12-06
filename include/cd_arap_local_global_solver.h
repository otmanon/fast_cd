#pragma once
#include "cd_sim_params.h"
#include "cd_arap_precomp.h"
#include "augment_with_linear_constraints.h"

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <igl/cat.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <json.hpp>
using namespace Eigen;
using namespace igl;
using namespace std;
using namespace nlohmann;
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

	SimplicialLDLT<SparseMatrix<double>>* ldlt_solver;

	cd_arap_local_global_solver() {};

	cd_arap_local_global_solver(SparseMatrix<double>& A, SparseMatrix<double>& Aeq, cd_arap_local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		SparseMatrix<double> H;
		augment_with_linear_constraints(A, Aeq, H);
		ldlt_solver = new SimplicialLDLT<SparseMatrix<double>>(H);
	}

	VectorXd solve(const VectorXd& z, const cd_sim_params& params, const cd_arap_dynamic_precomp& dp, const cd_arap_static_precomp& sp)
	{
		Eigen::VectorXd z_next = z, z_prev = z;
		if (p.to_convergence)
		{
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r = local_step(z_next, dp, sp);
				z_next = global_step(z_next, params, dp, sp, r);
				res = (z_next - z_prev).norm();
			} while (res > p.threshold);
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
		}
		return z_next;
	};

	VectorXd local_step(const VectorXd& z, const cd_arap_dynamic_precomp& dp, const cd_arap_static_precomp& sp)
	{
		VectorXd f = sp.K * z + dp.Kur + sp.Kx;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	
		
		Matrix3d F, rot;
		for (int i = 0; i < nt; i++)
		{
			F =  F_stack.block(3 * i, 0, 3, 3);// *sp.tet_vols(i);
		
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		return r;
	}

	VectorXd global_step(const VectorXd& z,const  cd_sim_params& params, const cd_arap_dynamic_precomp& dp, const cd_arap_static_precomp& sp, VectorXd& r)
	{
		VectorXd inertia_grad = -dp.My;
		VectorXd arap_grad = dp.Cur + sp.Cx - sp.VK.transpose() * r;

		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad + dp.f_ext;
		VectorXd Y;

		Eigen::VectorXd rhs = -igl::cat(1, g, dp.bc);
		VectorXd z_next = ldlt_solver->solve(rhs);
		
		return z_next.topRows(params.X.rows()*3);
	}

};