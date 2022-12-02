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
	Initializes an arap simulator that uses fast complementary dynamics.
	X -> mesh geometry
	T -> tet connectivity
	J -> control rig jacobian
	sim_dir -> directory to simulation cache, with a sim_params.json,
	solver_params.json, a clusters/folder a modes/ folder and a precomp/ folder

	If any of these are not found, will recompute everything from scratch assuming default parameters
	*/
	/*
	fast_cd_arap_sim(const MatrixXd& X, const MatrixXi& T, const SparseMatrix<double>& J, string& sim_dir)
	{
		//read sim_params.json, clusters and modes in this call
		params = fast_cd_sim_params(X, T, J, sim_dir);

		//read solver_params.json if not found, assume default
		cd_arap_local_global_solver_params& solver_params = cd_arap_local_global_solver_params(sim_dir + "/solver_params.json");

		if (!sp.read_from_cache(sim_dir + "/precomp/"))
		{
			printf("Could not read fast_cd_arap_sim matrix precomputations from cache %s, recomputing...\n", (sim_dir + "/precomp").c_str());
			sp = fast_cd_arap_static_precomp(params);
		}

		dp = fast_cd_arap_dynamic_precomp();
		sol = fast_cd_arap_local_global_solver(sp.BAB, sp.AeqB, solver_params);
	}*/

	/*
	initializes from cache dir. If not found, throws error. DOES NOT CHECK IF THE CACHE IS OTUDATED... do that somewhere else 
	*/
	fast_cd_arap_sim(std::string& cache_dir, fast_cd_sim_params& sim_params, cd_arap_local_global_solver_params& solver_params, bool read_cache, bool write_cache)
	{
		namespace fs = std::filesystem;

		this->params = sim_params;
		
		sp = fast_cd_arap_static_precomp();

		if (read_cache)
		{
			bool success = sp.read_from_cache(cache_dir);
			if (!success)
			{
				printf(" cache dir %s, is either corrupt or outdated. please construct fast_cd_arap_sim differently \n", cache_dir.c_str());
				printf(" Computing fast_cd_arap precomputations from scratch... \n", cache_dir.c_str());
				sp = fast_cd_arap_static_precomp(sim_params);
				if (write_cache)
					sp.write_to_cache(cache_dir);
				
			} 

		}
		else
		{
			printf(" Computing fast_cd_arap precomputations from scratch... \n", cache_dir.c_str());
			sp = fast_cd_arap_static_precomp(sim_params);
			if (write_cache)
				sp.write_to_cache(cache_dir);
			
		}

		dp = fast_cd_arap_dynamic_precomp();
		sol = fast_cd_arap_local_global_solver(sp.BAB, sp.AeqB, solver_params);
	}

	/*

	fast_cd_arap_sim(std::string& precomp_cache_dir, std::string& modes_cache_dir, std::string& clusters_cache_dir,
		cd_arap_local_global_solver_params& solver_params)
	{
		namespace fs = std::filesystem;

		params = fast_cd_sim_params();
		params.read_from_cache(modes_cache_dir, clusters_cache_dir);

		sp = fast_cd_arap_static_precomp();
		sp.read_from_cache(precomp_cache_dir);

		dp = fast_cd_arap_dynamic_precomp();
		sol = fast_cd_arap_local_global_solver(sp.BAB, sp.AeqB, solver_params);
	}*/

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
		assert(sp.AeqB.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp.BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z.rows() == sp.BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
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
		assert(params.Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp.BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z.rows() == sp.BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
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