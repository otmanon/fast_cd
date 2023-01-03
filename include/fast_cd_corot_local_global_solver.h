#pragma once
#include "fast_cd_corot_precomp.h"
#include "fast_cd_corot_sim_params.h"
#include "local_global_solver_params.h"
#include "augment_with_linear_constraints.h"

#include <igl/polar_svd3x3.h>
#include <igl/cat.h>
using namespace Eigen;
struct fast_cd_corot_local_global_solver
{
	local_global_solver_params sol_p;

	LLT<MatrixXd> llt_solver;
	LDLT<MatrixXd> ldlt_solver;

	fast_cd_corot_local_global_solver() {};
	/*
	Constructs arap local global solver object used to solve
	dynamics quickly for fast Complementary DYnamoics

	Inputs :
		A - m x m system matrix
		Aeq - c x m constraint rows that enforece Aeq z = b
				as linear equality constraints
		 sol_p - local_global_solver_params
	*/
	fast_cd_corot_local_global_solver(MatrixXd& A, MatrixXd& Aeq,
		local_global_solver_params& p)
	{
		this->sol_p = p;
		VectorXi bI;
		MatrixXd H;
		augment_with_linear_constraints(A, Aeq, H);

		if (Aeq.rows() > 0)
			ldlt_solver.compute(H);
		else
			llt_solver.compute(H);

	}

	//TODO this is the same code as alll the other local globals olvers... need 
	//method of abstracting the two.
	VectorXd solve(const VectorXd& z, const VectorXd& p, const fast_cd_corot_sim_params& params,
		const fast_cd_corot_dynamic_precomp& dp, const fast_cd_corot_static_precomp& sp)
	{
		Eigen::VectorXd z_next = z, z_prev = z;
		if (sol_p.to_convergence)
		{
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r, vg;
				local_step(z_next, p, dp, sp, r, vg);
				z_next = global_step(z_next, p, params, dp, sp, r, vg);
				res = (z_next - z_prev).norm();
			} while (res > sol_p.threshold);
		}
		else
		{
			int iter = 0;
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r, vg;
				local_step(z_next, p,  dp, sp, r, vg);
				z_next = global_step(z_next, p,  params, dp, sp, r, vg);
				res = (z_next - z_prev).norm();
				iter += 1;
				if (iter >= sol_p.max_iters)
					break;
			} while (res > sol_p.threshold);
		}
		return z_next;
	};


	//TOGO this is almost the same as the full space local global solver...
	void local_step(const VectorXd& z, const VectorXd& p, 
		const fast_cd_corot_dynamic_precomp& dp,
		const fast_cd_corot_static_precomp& sp,
		VectorXd& r, VectorXd& vol_g)
	{
		VectorXd f = sp.GKB * z + sp.GKJ*p;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		MatrixXd VG = MatrixXd::Zero(F_stack.rows(), F_stack.cols());

		Matrix3d F, rot, vol;

		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3).transpose();// *sp.tet_vols(i);

			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot.transpose();
			vol = (rot.transpose()* F - Matrix3d::Identity()).trace() * rot.transpose(); // trace(R'*F - eye(2))*R;
			VG.block(3 * i, 0, 3, 3) = vol;
		}
		r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		vol_g = Map<const VectorXd>(VG.data(), VG.rows() * VG.cols());
	}

	//TODO project z onto the constraint set if it isn't satisfied at first
	VectorXd global_step(const VectorXd& z,const VectorXd&p,
		const  fast_cd_corot_sim_params& params,
		const fast_cd_corot_dynamic_precomp& dp,
		const fast_cd_corot_static_precomp& sp,
		VectorXd& r, VectorXd& vol_g)
	{
		VectorXd inertia_grad = -dp.BMy;
		VectorXd arap_grad = sp.BCmuJ*p - sp.GVmuKB.transpose() * r;

		VectorXd vol_grad = sp.GVlamKB.transpose() * vol_g;
		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad
			+ vol_grad + dp.f_ext + sp.BAB * z;
		VectorXd Y;

		Eigen::VectorXd rhs = -igl::cat(1, g, dp.bc);

		VectorXd dz;
		if (params.AeqB.rows() > 0)
			dz = ldlt_solver.solve(rhs);
		else
			dz = llt_solver.solve(rhs);
		VectorXd z_next = z + dz.topRows(params.B.cols());

		////hould do a line search here:
		double alpha = 2.0;
		double E0 = fast_cd_energy(z, p, params, dp, sp);
		double E = E0;
		do
		{
			alpha /= 2.0;
			z_next = z + alpha * dz;
			E = fast_cd_energy(z_next, p, params, dp, sp);
		} while (E > E0 + 1e-8);

		std::cout << alpha << std::endl;
		return z_next;
	}


	double fast_cd_corot_energy(const VectorXd& z, const VectorXd& p,
		const  fast_cd_corot_sim_params& params,
		const fast_cd_corot_dynamic_precomp& dp,
		const fast_cd_corot_static_precomp& sp)
	{
		VectorXd f = sp.GKB * z + sp.GKJ * p;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Matrix3d F, rot;

		double arap_energy = 0;
		double vol_energy = 0;
		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3).transpose();// *sp.tet_vols(i);

			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot.transpose(); //rot needs to be transposed here to be correct

			Matrix3d temp = rot.transpose() * F - Matrix3d::Identity();
			vol_energy += 0.5 * sp.cluster_vols(i) * sp.cluster_lambda(i) *
				temp.trace() * temp.trace();
			// ;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		double arap_energy_efficient = (double)(0.5 * z.transpose() * (sp.BCmuB * z)) +
			z.transpose() * (sp.BCmuJ*p - sp.GVmuKB.transpose() * r);
		double total_energy = arap_energy_efficient +  vol_energy;
		return total_energy;
	}

	double fast_cd_kinetic_energy(const VectorXd& z, const VectorXd& p,
		const  fast_cd_corot_sim_params& params,
		const fast_cd_corot_dynamic_precomp& dp, const fast_cd_corot_static_precomp& sp)
	{
		double energy = params.invh2 * params.do_inertia *
			((0.5 * double(z.transpose() * sp.BMB * z)) - (double)(z.transpose() * dp.BMy));

		return energy;
	}

	double fast_cd_energy(const VectorXd& z,const VectorXd& p,
		const  fast_cd_corot_sim_params& params,
		const fast_cd_corot_dynamic_precomp& dp,
		const fast_cd_corot_static_precomp& sp)
	{
		double corot = fast_cd_corot_energy(z,p,  params, dp, sp);
		double inertia = fast_cd_kinetic_energy(z,p,  params, dp, sp);

		double ext = z.transpose() * dp.f_ext;
		double total_energy = corot + inertia + ext;
		return total_energy;
	};

};