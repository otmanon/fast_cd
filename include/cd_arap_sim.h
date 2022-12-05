#pragma once
#include "cd_sim_params.h"
#include "cd_sim_state.h"
#include "cd_arap_local_global_solver.h"
#include "cd_arap_precomp.h"



#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std; using namespace Eigen;



struct cd_arap_sim
{
public:

	cd_sim_params params;

	cd_arap_local_global_solver sol;

	cd_arap_static_precomp sp;

	cd_arap_dynamic_precomp dp;

	cd_arap_sim() {};
	cd_arap_sim(cd_sim_params& params, cd_arap_local_global_solver_params& solver_params) {
		this->params = params;

		this->sp = cd_arap_static_precomp(params);
		this->dp = cd_arap_dynamic_precomp();

		this->sol = cd_arap_local_global_solver(sp.A, params.Aeq, solver_params);

	};
	cd_arap_sim(cd_sim_params& params, cd_arap_local_global_solver_params& solver_params, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		this->params = params;

		this->sp = sp;
		this->dp = dp;
		this->sol = cd_arap_local_global_solver(sp.A, params.Aeq, solver_params);


	}

	virtual Eigen::VectorXd step(const VectorXd& z, const VectorXd& p, const cd_sim_state& state,const  VectorXd& f_ext,const  VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params.Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");

		dp.precomp(z, p, state, f_ext, bc, sp);
		VectorXd z_next = sol.solve(z, params, dp, sp);
		return z_next;
	}

	// no rig
	virtual Eigen::VectorXd step(const VectorXd& z, const cd_sim_state& state, const VectorXd& f_ext, const VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params.Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");

		dp.precomp(z, state, f_ext, bc, sp);
		VectorXd z_next = sol.solve(z, params, dp, sp);
		return z_next;
	}




	/*
	Updates the equality constraint matrix to be C, where C slices out known rows of the full space system 
	*/
	virtual void set_equality_constraint(SparseMatrix<double>& C)
	{
		this->params.Aeq = C;

		sp = cd_arap_static_precomp(this->params);
		cd_arap_local_global_solver_params old_params = sol.p;

		sol = cd_arap_local_global_solver(this->sp.A, this->params.Aeq, old_params);
	}



	/*
	Computes total kinetic energy
	*/
	virtual double kinetic_energy(VectorXd& z, VectorXd& p, cd_sim_state& state)
	{
		// Kinetic energy K = (x - y)^T M (x - y) , y = x_hist = 2x_curr - x_prev = 2(uc + ur + x0)_curr - (uc + ur + x0)_prev = 
		//															   2 uc_curr - uc_prev + 2ur_curr - ur_prev + x0
		//           K = (uc + ur + x0 - y)^T M (uc + ur + x0 - y)
		//           K = (uc + ur + x0 - ( 2 uc_curr - uc_prev + 2ur_curr - ur_prev + x0)) )^T M (uc + ur + x0 - ( 2 uc_curr - uc_prev + 2ur_curr - ur_prev + x0)))
		//			 K = (uc + ur  - ( 2 uc_curr - uc_prev + 2ur_curr - ur_prev)) )^T M (uc + ur + x0 - ( 2 uc_curr - uc_prev + 2ur_curr - ur_prev)))          // c = uc - 2uc_curr + uc_prev   // r = ur - 2ur_curr + ur_prev
		//			 K = (c + r)^TM(c+r) = c^T M c + r^T M r + 2 r^T M c;
		//			 K =  u'B'MBu +   (Jp - x0)'M(Jd - x0) + 2(Jd - x0)'M Bu                 //c = Bu  r=Jd - x0
		//			 K = u'B'MBu + d'J'MJd + x0'Mx0 - 2d'J'Mx0 + 2d'J'MBu - 2x0'MBu
		VectorXd u = z - 2 * state.z_curr + state.z_prev;
		VectorXd d = p - 2 * state.p_curr + state.p_prev;

		double ku = u.transpose() * sp.M * u;
		double kd = d.transpose() * sp.JMJ * d;
		double kx = sp.xMx(0);
		double kdx = -2 * d.transpose() * sp.JMx;
		double kdu = 2 * d.transpose() * sp.MJ.transpose() * u;
		double kxu = -2 * sp.Mx.transpose() * u;
		double total = params.invh2 * (ku + kd + kx + kdx + kdu + kxu);
		return total;
	}


	/*
	Computes kinetic energy, assuming zero rig displacement.
	*/
	virtual double kinetic_energy_complementary(VectorXd& z, cd_sim_state& state)
	{
		// Kinetic energy K = (x - y)^T M (x - y) , y = x_hist = 2x_curr - x_prev = 2(uc + ur + x0)_curr - (uc + ur + x0)_prev = 
		//																		   2 uc_curr - uc_prev + 2ur_curr - ur_prev + x0
		//                K = (uc + ur + x0 - y)^T M (uc + ur + x0 - y)
		//                K = (uc +  x0 - y)^T M (uc + x0 - y)  //assume no ur    y =  2 uc_curr - uc_prev  + x0
		//                K = (uc + x0 - (2uc_curr - uc_prev + x0))^T M (uc + x0 - (2uc_curr - uc_prev + x0)   
		//				  K = (uc - 2 uc_curr + uc_prev)^T M (uc - 2uc_curr +uc_prev)
		//                K = (z - 2z_curr + z_prev)^T B^T M B (z - 2z_curr + z_prev)
		VectorXd u = z - 2 * state.z_curr + state.z_prev;
		double k = params.invh2 * u.transpose() * sp.M * u;
		return k;
	}



	/*
	Projects z to full space
	*/
	virtual VectorXd full_space(VectorXd& z)
	{
		return z; //for a full space sim, z is already full space
	}
};