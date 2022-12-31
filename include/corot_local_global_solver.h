#pragma once
#include "sim_params.h"
#include "corot_precomp.h"
#include "augment_with_linear_constraints.h"
#include "local_global_solver_params.h"

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


struct corot_local_global_solver
{
	local_global_solver_params p;
	VectorXd z;

	SimplicialLDLT<SparseMatrix<double>>* ldlt_solver;

	corot_local_global_solver() {};

	corot_local_global_solver(SparseMatrix<double>& A, SparseMatrix<double>& Aeq, local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		SparseMatrix<double> H;
		augment_with_linear_constraints(A, Aeq, H);
		ldlt_solver = new SimplicialLDLT<SparseMatrix<double>>(H);
	}

	VectorXd solve(const VectorXd& z, const sim_params& params, const corot_dynamic_precomp& dp, const corot_static_precomp& sp)
	{
		Eigen::VectorXd z_next = z, z_prev = z;
		if (p.to_convergence)
		{
			double res;
			do
			{
				z_prev = z_next;
				VectorXd r, vg;
				local_step(z_next, dp, sp, r, vg);
				z_next = global_step(z_next, params, dp, sp, r, vg);
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
				VectorXd r, vg;
				local_step(z_next, dp, sp, r, vg);
				z_next = global_step(z_next, params, dp, sp, r, vg);
				res = (z_next - z_prev).norm();
				iter += 1;
				if (iter >= p.max_iters)
					break;
			} while (res > p.threshold);
		}
		return z_next;
	};

	void local_step(const VectorXd& z, const corot_dynamic_precomp& dp, const corot_static_precomp& sp,
		VectorXd& r, VectorXd& vol_g)
	{
		VectorXd f = sp.K * z + sp.Kx;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		MatrixXd VG = MatrixXd::Zero(F_stack.rows(), F_stack.cols());

		Matrix3d F, rot, vol;

		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);// *sp.tet_vols(i);

			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
			vol = (rot.transpose() * F - Matrix3d::Identity()).trace() * rot; // trace(R'*F - eye(2))*R;
			VG.block(3 * i, 0, 3, 3) = vol;
		}
		r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		
		vol_g = Map<const VectorXd>(VG.data(), VG.rows() * VG.cols());
	}

	VectorXd global_step(const VectorXd& z,
		const  sim_params& params, 
		const corot_dynamic_precomp& dp,
		const corot_static_precomp& sp, 
		VectorXd& r, VectorXd& vol_g)
	{
		VectorXd inertia_grad = -dp.My;
		VectorXd arap_grad = sp.Cmux - sp.VmuK.transpose() * r;

		VectorXd vol_grad = sp.VlamK.transpose() * vol_g;
		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad + vol_grad + dp.f_ext + sp.A * z;
		VectorXd Y;

		Eigen::VectorXd rhs = -igl::cat(1, g, dp.bc);
		VectorXd dz = ldlt_solver->solve(rhs);
		VectorXd z_next = z + dz.topRows(params.X.rows() * 3);


		//should do a line search here:
		double alpha = 2.0;
		double E0 = energy(z, params, dp, sp);
		double E = E0;
		do
		{
			alpha /= 2.0;
			z_next = z + alpha * dz;
			E = energy(z_next, params, dp, sp);
		} while (E > E0 + 1e-8);

		std::cout << alpha << std::endl;
		return z_next;
	}


	double corot_energy(const VectorXd& z,
		const  sim_params& params,
		const corot_dynamic_precomp& dp, const corot_static_precomp& sp)
	{
		VectorXd f = sp.K * z + sp.Kx;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Matrix3d F, rot;

		double arap_energy = 0;
		double vol_energy = 0;
		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);// *sp.tet_vols(i);

			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;// .transpose();

			Matrix3d temp = rot.transpose() * F - Matrix3d::Identity();
			vol_energy += sp.tet_vols(i) * params.lambda(i) *
				temp.trace() * temp.trace();
			// ;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());

		double arap_energy_efficient = 0.5 * (double)(z.transpose() * (sp.Cmu * z)) +
			z.transpose() * (sp.Cmux - sp.VmuK.transpose() * r);


		double total_energy = arap_energy_efficient + vol_energy;
		return total_energy;
	}

	double kinetic_energy(const VectorXd& z,
		const  sim_params& params,
		const corot_dynamic_precomp& dp, const corot_static_precomp& sp)
	{
		double energy = params.invh2 * params.do_inertia * 
		(0.5*double(z.transpose() *sp.M*z) 
			- z.transpose() * dp.My);
		return energy;
	}


	
	double energy(const VectorXd& z,
		const  sim_params& params, 
		const corot_dynamic_precomp& dp, const corot_static_precomp& sp)
	{
		double corot= corot_energy(z, params, dp, sp);
		double inertia = kinetic_energy(z, params, dp, sp);

		double ext = z.transpose() * dp.f_ext;
		double total_energy = corot + inertia + ext;
		return total_energy;
	};

};