#pragma once
#include "cd_sim_params.h"
#include "read_clusters_from_cache.h"
#include "read_modes_from_cache.h"

struct fast_cd_sim_params : cd_sim_params
{
	MatrixXd B;
	VectorXi labels;
	fast_cd_sim_params() {};
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T, 
		const MatrixXd& B, const VectorXi&  l,
		const SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia) : cd_sim_params(X, T, J, mu, lambda, h, do_inertia)
	{
		this->B = B;
		this->labels = l;

		//usually no equality constraint.
		this->Aeq.resize(0, X.rows() * 3);
	}


	/*
	Constructor that attempts to read skinning modes and clusters from cache. 
	If skinning modes or clusters can't be found, print and do nothing
	*/
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const string& mode_cache_dir, const string& clusters_cache_dir,
		const SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia) : cd_sim_params(X, T, J, mu, lambda, h, do_inertia)
	{
		read_from_cache(mode_cache_dir, clusters_cache_dir);
		//usually no equality constraint.
		this->Aeq.resize(0, X.rows() * 3);
	}

	bool read_from_cache(const string& mode_cache_dir, const string& clusters_cache_dir)
	{
		VectorXd L;
		bool read_modes = read_skinning_modes_from_cache(mode_cache_dir, this->B, L);
		bool read_clusters = read_clusters_from_cache(clusters_cache_dir, this->labels);

		if (!read_modes)
		{
			printf("Could not read skinning modes from cache directory %s \n", mode_cache_dir.c_str());
		}
		if (!read_clusters)
			printf("Could not read clusters  from cache directory %s \n", clusters_cache_dir.c_str());

		return read_modes && read_clusters;
	}

};