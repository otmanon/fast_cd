#pragma once
#include "cd_arap_sim.h"
#include "fast_cd_sim_params.h"
#include "fast_cd_arap_precomp.h"
#include "fast_cd_arap_local_global_solver.h"
#include "read_fast_cd_sim_static_precompute.h"
#include "write_fast_cd_sim_static_precomputation.h"


struct fast_cd_arap_sim : cd_arap_sim
{
public:

	fast_cd_sim_params params;

	fast_cd_arap_local_global_solver sol;

	fast_cd_arap_static_precomp sp;

	fast_cd_arap_dynamic_precomp dp;

	fast_cd_arap_sim() {};

	/*
	initializes from cache dir. If not found, throws error. DOES NOT CHECK IF THE CACHE IS OTUDATED... do that somewhere else 
	*/
	fast_cd_arap_sim(std::string& cache_dir, cd_arap_local_global_solver_params& solver_params)
	{
		namespace fs = std::filesystem;

		MatrixXd B; VectorXd L; VectorXi l;
		bool well_read = read_fast_cd_sim_static_precomputation(cache_dir, B, L, l, sp.BCB, sp.BMB, sp.BAB,
			sp.AeqB, sp.GmKB, sp.GmKJ, sp.GmKx, sp.G1VKB, sp.BMJ, sp.BMx, sp.BCJ, sp.BCx);

		if (!well_read)
			printf(" cache dir %s, is either corrupt or outdated. please construct fast_cd_arap_sim differently \n", cache_dir);

		params = fast_cd_sim_params();
		params.B = B;
		params.labels = l;

		dp = fast_cd_arap_dynamic_precomp();
		sol = fast_cd_arap_local_global_solver(sp.BAB, sp.AeqB, solver_params);
	}

	fast_cd_arap_sim(fast_cd_sim_params& sim_params, cd_arap_local_global_solver_params& solver_params) {
		sp =  fast_cd_arap_static_precomp(sim_params);
		dp =  fast_cd_arap_dynamic_precomp();
		params = sim_params;

		sol =  fast_cd_arap_local_global_solver(sp.BAB, sp.AeqB, solver_params);
	};
	fast_cd_arap_sim(fast_cd_sim_params& params, fast_cd_arap_local_global_solver& solver, fast_cd_arap_dynamic_precomp& dp, fast_cd_arap_static_precomp& sp) : cd_arap_sim()
	{
		this->params = params;
		this->sol = solver;
		this->sp = sp;
		this->dp = dp;
	}

	Eigen::VectorXd step(const VectorXd& z,  const VectorXd& p,  const cd_sim_state& state, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		////needs to be updated every timestep
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp->BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z_ext.rows() == sp->BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
		VectorXd z_next = z;
		dp.precomp(z, p, state, f_ext, bc, sp);
		z_next = sol.solve(z_next, params, dp, sp);
		return z_next;
	}

	Eigen::VectorXd step_test(const VectorXd& z, const VectorXd& p, const cd_sim_state& state, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		////needs to be updated every timestep
		std::cout << "went through this function flawlessly" << std::endl;
		return VectorXd();
	}

	Eigen::VectorXd step(const VectorXd& z, const VectorXd& p, const VectorXd& z_curr, const VectorXd& z_prev,
		const VectorXd& p_curr, const VectorXd& p_prev, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		////needs to be updated every timestep
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp->BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z_ext.rows() == sp->BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
		VectorXd z_next = z;
		cd_sim_state state(z_curr, z_prev, p_curr, p_prev);
		dp.precomp(z, p, state, f_ext, bc, sp);
		z_next = sol.solve(z_next, params, dp, sp);
		return z_next;
	}


	/*Saves all precomputed matrices required for simulation in cache_dir*/
	bool save(std::string& cache_dir)
	{
		VectorXd L; // L is useless at this point.
		bool well_saved = write_fast_cd_sim_static_precomputation(cache_dir, params.B, L, params.labels, sp.BCB, sp.BMB, sp.BAB,
			sp.AeqB, sp.GmKB, sp.GmKJ, sp.GmKx, sp.G1VKB, sp.BMJ, sp.BMx, sp.BCJ, sp.BCx);
	
		return well_saved;
	}


};