#pragma once
#include "cd_arap_sim.h"
#include "fast_cd_arap_sim_params.h"
#include "fast_cd_arap_precomp.h"
#include "fast_cd_arap_local_global_solver.h"
#include "read_fast_cd_sim_static_precompute.h"
#include "write_fast_cd_sim_static_precomputation.h"


struct fast_cd_arap_sim
{
public:

	
	fast_cd_arap_sim_params* params;

	fast_cd_arap_local_global_solver* sol;

	fast_cd_arap_static_precomp* sp;

	fast_cd_arap_dynamic_precomp* dp;

	fast_cd_arap_sim() {};

	/*
	initializes from cache dir. If not found,
	throws error. 
	*/
	fast_cd_arap_sim(std::string& cache_dir, fast_cd_arap_sim_params& sim_params, 
		local_global_solver_params& solver_params, bool read_cache, 
		bool write_cache)
	{
		namespace fs = std::filesystem;
		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp();
		if (read_cache)
		{
			bool success = fcd_sp->read_from_cache(cache_dir);
			if (!success)
			{
				printf(" cache dir %s, is either corrupt or outdated. please construct fast_cd_arap_sim differently \n", cache_dir.c_str());
				printf(" Computing fast_cd_arap precomputations from scratch... \n", cache_dir.c_str());
				fcd_sp = new fast_cd_arap_static_precomp(sim_params);
				/*if (write_cache)
					fcd_sp->write_to_cache(cache_dir);*/
			} 
		}
		else
		{
			printf(" Computing fast_cd_arap precomputations from scratch... \n", cache_dir.c_str());
			fcd_sp = new fast_cd_arap_static_precomp(sim_params);
			if (write_cache)
				fcd_sp->write_to_cache(cache_dir);
			
		}

		params = &sim_params;
		sp = fcd_sp;
		dp = new  fast_cd_arap_dynamic_precomp();
		sol =  new fast_cd_arap_local_global_solver(((fast_cd_arap_static_precomp*)sp)->BAB, 
			((fast_cd_arap_static_precomp*)sp)->AeqB, solver_params);

	}

	fast_cd_arap_sim(fast_cd_arap_sim_params& sim_params, local_global_solver_params& solver_params) 
	{
		params = &sim_params;
		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp(sim_params);

		sp = fcd_sp; // fast_cd_arap_static_precomp(sim_params);
		dp =  new fast_cd_arap_dynamic_precomp();

		sol = new fast_cd_arap_local_global_solver(fcd_sp->BAB, fcd_sp->AeqB, solver_params);
	};

	fast_cd_arap_sim(fast_cd_arap_sim_params& sim_params)
	{
		params = &sim_params;
		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp(sim_params);
		local_global_solver_params solver_params =  local_global_solver_params(false, 10, 1e-6);
		sp = fcd_sp; // fast_cd_arap_static_precomp(sim_params);
		dp = new fast_cd_arap_dynamic_precomp();

		sol = new fast_cd_arap_local_global_solver(fcd_sp->BAB, fcd_sp->AeqB, solver_params);
	};


	fast_cd_arap_sim(fast_cd_arap_sim_params& params, fast_cd_arap_local_global_solver& solver, fast_cd_arap_dynamic_precomp& dp, fast_cd_arap_static_precomp& sp)
	{
		this->params = &params;
		this->sol = &solver;
		this->sp = &sp;
		this->dp = &dp;
	}

	/*
	Advances the pre-configured simulation one step
	Inputs : 
		z :  m x 1 current guess for z ("maybe shouldn't expose this")
		sol_p :  sol_p x 1 flattened rig parameters following writeup column order flattening converntion
		state : sim_cd_state that contains info like z_curr, z_prev, p_curr and p_prev
		f_ext : used to specify excternal forces like gravity.
		bc :	rhs of equality constraint if some are configured in system
				(should match in rows with sim.params.Aeq)
	Outpus :
		z_next : m x 1 next timestep degrees of freedom
	*/
   virtual Eigen::VectorXd step(const VectorXd& z,  const VectorXd& p,  const cd_sim_state& state, const  VectorXd& f_ext, const  VectorXd& bc)
	{
	   fast_cd_arap_static_precomp* sp = (fast_cd_arap_static_precomp*)this->sp;
		////needs to be updated every timestep
		assert(sp->AeqB.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp->BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z.rows() == sp->BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
		VectorXd z_next = z;
		((fast_cd_arap_dynamic_precomp*)dp)->precomp(z, p, state, f_ext, bc, *sp);
		z_next = ((fast_cd_arap_local_global_solver*)sol)->solve(z_next, *((fast_cd_arap_sim_params*)params), *((fast_cd_arap_dynamic_precomp*)dp), *sp);
		return z_next;
	}

	

   /*
   Advances the pre-configured simulation one step
   Inputs :
	   z :  m x 1 current guess for z ("maybe shouldn't expose this")
	   sol_p :  sol_p x 1 flattened rig parameters following writeup column order flattening converntion
	   z_curr : m x 1 current d.o.f.s
	   z_prev : m x 1 previos d.o.f.s
	   p_curr : sol_p x 1 current rig parameters
	   p_prev : sol_p x 1 previous rig parameters
	   f_ext : used to specify excternal forces like gravity.
	   bc :	rhs of equality constraint if some are configured in system
			   (should match in rows with sim.params.Aeq)
   Outpus :
	   z_next : m x 1 next timestep degrees of freedom
   */
	virtual Eigen::VectorXd step(const VectorXd& z, const VectorXd& p, const VectorXd& z_curr, const VectorXd& z_prev,
		const VectorXd& p_curr, const VectorXd& p_prev, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		////needs to be updated every timestep
		fast_cd_arap_static_precomp* sp = (fast_cd_arap_static_precomp*)this->sp;
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp->BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z.rows() == sp->BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
		VectorXd z_next = z;
		cd_sim_state state(z_curr, z_prev, p_curr, p_prev);
		dp->precomp(z, p, state, f_ext, bc, *sp);
		z_next = sol->solve(z_next, *((fast_cd_arap_sim_params*)params), *((fast_cd_arap_dynamic_precomp*)dp), *sp);
		return z_next;
	}


	///*Saves all precomputed matrices required for simulation in cache_dir*/
	//bool save(std::string& cache_dir)
	//{
	//	VectorXd L; // L is useless at this point.
	//	bool well_saved = write_fast_cd_sim_static_precomputation(cache_dir, params.B, L, params.labels, sp.BCB, sp.BMB, sp.BAB,
	//		sp.AeqB, sp.GmKB, sp.GmKJ, sp.GmKx, sp.G1VKB, sp.BMJ, sp.BMx, sp.BCJ, sp.BCx);
	//
	//	return well_saved;
	//}


	/*
	Computes kinetic energy, assuming zero rig displacement.
	*/
	double kinetic_energy_complementary(VectorXd& z,  cd_sim_state& state)
	{
		// Kinetic energy K = (x - y)^T M (x - y) , y = x_hist = 2x_curr - x_prev = 2(uc + ur + x0)_curr - (uc + ur + x0)_prev = 
		//																		   2 uc_curr - uc_prev + 2ur_curr - ur_prev + x0
		//                K = (uc + ur + x0 - y)^T M (uc + ur + x0 - y)
		//                K = (uc +  x0 - y)^T M (uc + x0 - y)  //assume no ur    y =  2 uc_curr - uc_prev  + x0
		//                K = (uc + x0 - (2uc_curr - uc_prev + x0))^T M (uc + x0 - (2uc_curr - uc_prev + x0)   
		//				  K = (uc - 2 uc_curr + uc_prev)^T M (uc - 2uc_curr +uc_prev)
		//                K = (z - 2z_curr + z_prev)^T B^T M B (z - 2z_curr + z_prev)
		fast_cd_arap_static_precomp* sp = (fast_cd_arap_static_precomp*)this->sp;
		VectorXd u = z - 2 * state.z_curr + state.z_prev;
		double k = params->invh2 * u.transpose() * sp->BMB * u;
		return k;
	}

	/*
	Computes total kinetic energy
	*/
	double kinetic_energy(VectorXd& z, VectorXd& p, cd_sim_state& state)
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
		fast_cd_arap_static_precomp* sp = (fast_cd_arap_static_precomp*)this->sp;
		double ku = u.transpose() * sp->BMB * u;
		double kd = d.transpose() * sp->JMJ * d;
		double kx = sp->xMx(0);
		double kdx = -2 * d.transpose() * sp->JMx;
		double kdu = 2 * d.transpose() * sp->BMJ.transpose() * u;
		double kxu = -2 * sp->BMx.transpose() * u;
		double total = params->invh2 *( ku + kd + kx + kdx + kdu + kxu);
		return total;
	}


	/*
	Updates the equality constraints to new ones as given by matrix Aeq
	*/
	virtual void set_equality_constraint(SparseMatrix<double>& Aeq)
	{
		params->set_equality_constraint(Aeq);
		sol = new fast_cd_arap_local_global_solver(
			sp->A, params->Aeq,
			sol->p);
	}


	///*
	//Copmutes total arap energy
	//*/
	//double arap_energy(VectorXd& z, VectorXd& sol_p, cd_sim_state& state)
	//{
	//	//rotation fitting, then evaluate energy

	//}
	//
	///*
	//Copmutes the total energy (elastic + inertia) for query state z, rig parameters sol_p, and current state variables in state
	//*/
	//double energy(VectorXd& z, VectorXd& sol_p, cd_sim_state& state)
	//{
	//	double arap = arap_energy(z, sol_p, state);
	//	double kinetic = kinetic_energy(z, sol_p, state);
	//	double total = kinetic + arap;
	//	return total;
	//}

	VectorXd full_space(VectorXd& z)
	{
		fast_cd_arap_sim_params* params = (fast_cd_arap_sim_params*)this->params;
		VectorXd u = params->B * z;
		return u;
	}


	fast_cd_arap_sim_params* parameters()
	{
		return (fast_cd_arap_sim_params*) this->params;
	}


};