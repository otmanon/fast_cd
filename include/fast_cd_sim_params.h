#pragma once
#include "cd_sim_params.h"
#include "read_clusters_from_cache.h"
#include "read_modes_from_cache.h"
#include "lbs_jacobian.h"

using namespace Eigen;
using namespace std;
void default_sim_params(double& mu, double& lambda, double& h, bool& do_inertia)
{
	mu = 10;
	lambda = 0;
	h = 1e-2;
	do_inertia = true;
}

struct fast_cd_sim_params : cd_sim_params
{
	MatrixXd B;
	VectorXi labels;

	string mode_type; // skinning or displacement.
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


	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T, const SparseMatrix<double>& J, const MatrixXd& B, const VectorXi& labels)
	{
		this->X = X;
		this->T = T;
		this->J = J;

		this->B = B;
		this->labels = labels;
		double mu; double lambda;
		default_sim_params(mu, lambda, this->h, this->do_inertia);
		this->mu = mu * VectorXd(T.rows(), 1);
		this->lambda = lambda * VectorXd(T.rows(), 1);
		this->invh2 = 1.0 / (h * h);
	};

	/*
	Constructor that attempts to read skinning modes and clusters from cache. 
	If skinning modes or clusters can't be found, print and do nothing
	*/
	/*
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const string& mode_cache_dir, const string& clusters_cache_dir,
		const SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia) : cd_sim_params(X, T, J, mu, lambda, h, do_inertia)
	{
		read_from_cache(mode_cache_dir, clusters_cache_dir);
		//usually no equality constraint.
		this->Aeq.resize(0, X.rows() * 3);
	};
	*/

	
	/*
	Constructor that builds a fast cd arap sim from scratch by loading up simulation parameters from cache
	*/
	/*
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T, const SparseMatrix<double>& J, const string& params_dir)
	{
		this->X = X;
		this->T = T;
		this->J = J;

		if (!this->read_from_cache(params_dir))
		{
			double mu; double lambda;
			default_sim_params(mu, lambda, this->h, this->do_inertia);
			this->mu = mu * VectorXd(T.rows(), 1);
			this->lambda = lambda * VectorXd(T.rows(), 1);
			this->invh2 = 1.0 / (h * h);
		}


	}
	*/
	/*
	Reads simulation parameters from scratch
	
	bool read_from_cache(const string& params_dir)
	{

		bool success = true;

		//read params
		double mu; double lambda;
		success = success && read_sim_params_from_cache(params_dir + "/sim_params.json", mu, lambda, h, do_inertia);
		this->mu = mu * VectorXd(T.rows(), 1);
		this->lambda = lambda * VectorXd(T.rows(), 1);
		this->invh2 = 1.0 / (h * h);

		// read modes
		if (mode_type == "skinning")
		{
			MatrixXd W;
			VectorXd L;
			SparseMatrix<double> Js;
			success = read_skinning_modes_from_cache(params_dir + "/modes/", W, L);
			lbs_jacobian(X, W, Js);
			B = Js.toDense();
		}
		else if (mode_type == "displacement")
		{
			VectorXd L;
			success = read_displacement_modes_from_cache(params_dir + "/modes/", B, L);
		}
		//read clusters
		success = success && read_clusters_from_cache(params_dir + "/clusters/", this->labels);
		return success;
	}
	*/

	/*
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
	}*/

};