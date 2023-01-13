#pragma once
#include "lbs_jacobian.h"
#include "create_directory_by_force.h"
#include "compute_clusters_igl.h"

#include "get_modes.h"

#include <Eigen/Core>
#include <igl/readDMAT.h>

#include <filesystem>
using namespace Eigen;
using namespace std;
struct fast_cd_subspace_parameters
{
	int num_modes;
	int num_clusters;

	int num_clustering_features;

	string mode_type;
	string subspace_constraint_type;
	bool split_components;

	//for development resons

	bool debug;
	string output_dir;
	//default constructor
	fast_cd_subspace_parameters() {}
	fast_cd_subspace_parameters(int num_modes, string& subspace_constraint_type, string& mode_type, int num_clusters, int num_clustering_features, bool split_components, bool debug=false,
		string output_dir="") {
		this->num_modes = num_modes;
		this->num_clusters = num_clusters;
		this->num_clustering_features = num_clustering_features;
		this->mode_type = mode_type;
		this->split_components = split_components;
		this->subspace_constraint_type = subspace_constraint_type;
		this->debug = debug;
		this->output_dir = output_dir;


		if (debug)
		{
			if (output_dir == "")
				printf("output_dir is blank, have no where to write debug info! \n");
			else
				create_directory_by_force(output_dir);
		}
	}

};
/*
This class is in charge of constructing a subspace

*/
struct  fast_cd_subspace
{

	fast_cd_subspace_parameters params;

	MatrixXd B;    // final subspace
	MatrixXd W;    //skinning weights (if applicable)
	VectorXd L;   // eigenvalues

	VectorXi l;   // clusters

	SparseMatrix<double> Aeq; //custom constraint. empty if not using custom constraint c x 3V

	fast_cd_subspace() {};
	/*
	Initialized and configures subspace... does NOT build it yet.
	Inputs:
	num_modes - (int) number of modes in subspace
	subspace_constraint_type - (string) either "none", or "cd" or "cd_momentum_leak"
	mode_type - (string) either "skinning" or "displacement"
	num_clusters - (int) number of clusters
	num_clustering_features - (int) number of clustering features used
	split_components - (bool) whether to split the components of the clusters 

	Optional:
	debug - (bool) whether to save and store debug info 
	output_dir - (string) output directory where we will write debug info
	*/
	fast_cd_subspace(int num_modes, string subspace_constraint_type, string& mode_type, int num_clusters, int num_clustering_features, bool split_components, bool debug = false,
		string output_dir = "")
	{
		params = fast_cd_subspace_parameters(num_modes, subspace_constraint_type, mode_type, num_clusters, num_clustering_features, split_components, debug, output_dir);
	}

	fast_cd_subspace(fast_cd_subspace_parameters& p)
	{
		params = p;
	};
	//

	/*
	Computes modes + clusters from scratch
	Inputs :
	V -> |n|x3 geometry
	T -> |T|x4 tet indices
	J -> |c|x3n null space/linear orthogonality constraint
	*/
	void init(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& J)
	{
		printf("Computing %i modes from scratch, this will cause Matlab to open... \n", params.num_modes);
		get_modes(V, T, J, params.mode_type, params.num_modes, B, L, W, params.debug, params.output_dir);
		printf("Computing %i clusters from scratch... \n", params.num_clusters);
		MatrixXd C;// cluster centers in weight space... just throw this out//
		if (params.mode_type == "displacement")
		{
			compute_clusters_displacement_features(T, B, L, params.num_clusters, params.num_clustering_features, l, C, false);
		}
		else if (params.mode_type == "skinning")
		{
			compute_clusters_weight_features(T, W, L, params.num_clusters, params.num_clustering_features, l, C, false);
		}
	}

	/*
	Computes modes + clusters from scratch, with control over when to read/write from cache
	Inputs:
	V -> |n|x3 geometry
	T -> |T|x4 tet indices
	J -> |c|x3n null space/linear orthogonality constraint
	
	read_cache -> whether or not to attempt to read modes and clusters from cache.
	write_cache -> whether or not to write modes and clusters to cache.
	modes_cache_dir -> directory where mode cache is
	(Optional)
	clusters_cache_dir -> directory where clusters cache is
	recompute_modes_if_not_found -> whether or not to recompute modes from scratch if not found in cache (default true)
	recompute_clusters_if_not_found -> whether or not to recompute clusters from scratich if not found in cache (default true)
	*/
	void init_with_cache(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& J,  bool read_cache, bool write_cache,
		string& modes_cache_dir, string& clusters_cache_dir, bool recompute_modes_if_not_found = true,
		bool recompute_clusters_if_not_found = true)
	{
		SparseMatrix<double> Aeq = SparseMatrix<double>(0, V.rows()*V.cols());
		if (params.subspace_constraint_type == "cd")
		{
			SparseMatrix<double>  M;
			igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
			Aeq = (J.transpose() * igl::repdiag(M, 3));
		}
		else if (params.subspace_constraint_type == "cd_momentum_leak")
		{
			SparseMatrix<double> D, M;
			momentum_leaking_matrix(V, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
			igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
			Aeq = (J.transpose() * igl::repdiag(M, 3)*igl::repdiag(D, 3));
		}
		else if (params.subspace_constraint_type == "custom")
		{
			if (this->Aeq.rows() == 0 || this->Aeq.cols() != V.rows()*3)
			{
				printf("Aeq matrix does not have the right amount of rows or columns to be valid. Please set the constraint \n \
				  by first calling fast_cd_subspace::set_custom_subspace_constraint(Aeq) before this call \n");
			}
			Aeq = this->Aeq;
		}
		// load subspace from cache... if modes can't be found, they are recomputed. if clusters can't be found, they are recomputed
		if (read_cache)
		{
			bool found_subspace_cache = this->read_from_cache_recompute(V, T, Aeq, modes_cache_dir, clusters_cache_dir,
				recompute_modes_if_not_found, recompute_clusters_if_not_found, write_cache);
			//if (write_cache && !found_subspace_cache)
			//	write_to_cache(modes_cache_dir, clusters_cache_dir);
		}
		else
		{
			this->init(V, T, Aeq);
			if (write_cache)
				this->write_to_cache(modes_cache_dir, clusters_cache_dir);
		}
	}

	
	void set_custom_subspace_constraint(SparseMatrix<double>& Aeq)
	{
		this->Aeq = Aeq;
	}

	/*
	Reads modes, clusters from cache directories. If they are not found, recompute each one of these as specified by
	recomp_modes_if_not_found//
	*/
	bool read_from_cache_recompute(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& J,
		std::string& modes_dir, std::string& clusters_dir, bool recomp_modes_if_not_found,
		bool recomp_clusters_if_not_found, bool write_cache)
	{
		bool success = true;
		success = success && read_modes_from_cache(modes_dir);
		
		if (success && params.mode_type == "skinning")
		{
			SparseMatrix<double> Bs;
			lbs_jacobian(V, W, Bs);
			B = Bs.toDense();
		}
		if (!success && recomp_modes_if_not_found)
		{
			printf("Computing %i modes from scratch, this will cause Matlab to open... \n", params.num_modes);
			get_modes(V, T, J, params.mode_type, params.num_modes, B, L, W, params.debug, params.output_dir);
			if (write_cache)
				write_modes_to_cache(modes_dir);
		}
		success = success && read_clusters_from_cache(clusters_dir);
		if (!success && recomp_clusters_if_not_found)
		{
			printf("Computing %i clusters from scratch... \n", params.num_clusters);
			MatrixXd C;// cluster centers in weight space... just throw this out
			if (params.mode_type == "displacement")
			{
				compute_clusters_displacement_features(T, B, L, params.num_clusters, params.num_clustering_features, l, C, params.split_components);
			}
			else if (params.mode_type == "skinning")
			{
				//compute_clusters_displacement_features(T, B, L, params.num_clusters, params.num_clustering_features, l, C, true);
				compute_clusters_weight_features(T, W, L, params.num_clusters, params.num_clustering_features, l, C, params.split_components);
			}
			if (write_cache)
				write_clusters_to_cache(clusters_dir);
		}
		return success;
	}

	/*
	Writes modes and clusters to cache directories
	modes_dir -> (string) where to save modes directory (both B.DMAT and L.DMAT, modes + frequencies)
	clusters_dir -> (string) where to save clusters directory (labels.DMAT)
	*/
	bool write_to_cache(std::string& modes_dir, std::string& clusters_dir)
	{
		bool success = true;
		success = success && write_modes_to_cache(modes_dir);
		success = success && write_clusters_to_cache(clusters_dir);
		return success;
	}

	bool write_modes_to_cache(std::string& modes_dir)
	{
		namespace fs = std::filesystem;
		if (!fs::exists(fs::path(modes_dir)))
		{
			fs::create_directories(fs::path(modes_dir));
		}
		printf("Writing modes to cache dir %s ... \n", modes_dir.c_str());
		bool success = true;
		if (params.mode_type == "skinning")
			success = success && igl::writeDMAT(modes_dir + "W.DMAT", W);
		else if (params.mode_type == "displacement")
			success = success && igl::writeDMAT(modes_dir + "B.DMAT", B);
		success = success && igl::writeDMAT(modes_dir + "L.DMAT", L);
		if (!success) printf("Could not save modes cache!");
		return success;
	}

	bool write_clusters_to_cache(std::string& clusters_dir)
	{
		namespace fs = std::filesystem;
		if (!fs::exists(fs::path(clusters_dir)))
		{
			fs::create_directories(fs::path(clusters_dir));
		}
		printf("Writing clusters to cache dir %s ... \n", clusters_dir.c_str());
		bool success = true;
		success = success && igl::writeDMAT(clusters_dir + "labels.DMAT", l);
		if (!success) printf("Could not save clusters cache!");
		return success;
	}


	/*
	Reads modes and clusters from cache directories
	Inputs:
	  modes_dir - (string) directory where B.DMAT/W.DMAT and L.DMAT is stored
	  clusters_dur - (string) directory where cluster labels.DMAT is stored
	*/
	bool read_from_cache(std::string& modes_dir, std::string& clusters_dir)
	{
		bool success = true;
		success = success && read_modes_from_cache(modes_dir);
		success = success && read_clusters_from_cache(clusters_dir);
		return success;
	}


	/*
	Read clusters from cache directories
	Inputs:
	  clusters_dir - (string) directory where cluster labels.DMAT is stored
	*/
	bool read_clusters_from_cache(std::string& clusters_dir)
	{
		bool success = true;
		success = success && igl::readDMAT(clusters_dir + "/labels.DMAT", l);
		bool correct_num_clusters = true;
		//allow for some wiggle room in cluster spec
		if (success)
			bool correct_num_clusters =  l.maxCoeff()+1 < params.num_clusters + 10 && l.maxCoeff()+1 > params.num_clusters - 10;
		
		success = success && correct_num_clusters;
		if (success)
		{
			printf("Successfully Read %s Clusters from cache dir %s! \n", params.mode_type.c_str(), clusters_dir.c_str());
		}
		if (!success)
		{
			if (!correct_num_clusters) printf("Num clusters in cache (%i) does not match num clusters requested (%i)", l.maxCoeff() + 1, params.num_clusters);
			printf("Could not read %s Clusters from cache dir %s! \n", params.mode_type.c_str(), clusters_dir.c_str());
		}
		return success;
	}
	
	/*
	Read clusters from cache directories
	Inputs:
	  modes_dur - (string) directory where cluster B.DMAT/W.DMAT and L.DMAT is stored
	*/
	bool read_modes_from_cache(std::string& modes_dir)
	{
		bool success = true;
		printf("Reading %s Modes from cache dir %s ... \n", params.mode_type.c_str(), modes_dir.c_str());
		bool correct_num_modes = true;
		int found_num_modes = 0;
		if (params.mode_type == "skinning")
		{
			success = success && igl::readDMAT(modes_dir + "/W.DMAT", W);
			success = success && igl::readDMAT(modes_dir + "/L.DMAT", L);
		
			if (success)
				found_num_modes = W.cols();
			
		}
		else if (params.mode_type == "displacement")
		{
			success = success && igl::readDMAT(modes_dir + "/B.DMAT", B);
			success = success && igl::readDMAT(modes_dir + "/L.DMAT", L);
			if (success)
				found_num_modes = B.cols();
		}
		if (success)
		{
			correct_num_modes = found_num_modes >= params.num_modes;
		}
		success = success && correct_num_modes;
		if (success)
		{
			if (params.mode_type == "skinning")
				W = W.leftCols(params.num_modes).eval();
			else if (params.mode_type == "displacement")
				B = B.leftCols(params.num_modes).eval();
			L = L.topRows(params.num_modes);
			

			printf("Successfully Read %s Modes from cache dir %s! \n", params.mode_type.c_str(), modes_dir.c_str());
		}
		if (!success)
		{
			if (!correct_num_modes) printf("Num modes in cache (%i) does not match num modes requested (%i)", found_num_modes, params.num_modes);
			printf("Could not read %s Modes from cache dir %s! \n", params.mode_type.c_str(), modes_dir.c_str());
		}
		return success;
	}
};