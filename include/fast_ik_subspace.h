#pragma once

#include "create_directory_by_force.h"
#include "compute_clusters_igl.h"

#include <Eigen/Core>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>

#include <filesystem>

using namespace Eigen;
using namespace std;

/*
This class is in charge of constructing a subspace

*/
struct  fast_ik_subspace
{

	int num_clusters;
	VectorXi l;   // clusters

	fast_ik_subspace() {};
	/*
	Initialized and configures fast_ik subspace, which invovles only computing clusters
	does NOT build it yet.
	Inputs:
	num_clusters - (int) number of clusters

	Optional:
	debug - (bool) whether to save and store debug info
	output_dir - (string) output directory where we will write debug info
	*/
	fast_ik_subspace( int num_clusters, bool debug = false,
		string output_dir = "")
	{
		this->num_clusters = num_clusters;
	}

	
	//

	/*
	Computes  clusters from scratch with input weight features at each vertex
	Inputs :
	V -> |n|x3 geometry
	T -> |T|x4 tet indices
	W -> |n| x m Clustering features
	*/
	void init(MatrixXd& V, MatrixXi& T, MatrixXd& W)
	{
		printf("Computing %i clusters from scratch... \n", num_clusters);
		MatrixXd C;// cluster centers in weight space... just throw this out//
		int m = W.cols();
		VectorXd L = VectorXd::Ones(m); //importance of each feature

		compute_clusters_weight_features( T, W, L, num_clusters, m, l, C, false);
	
	}

	/*
	Computes clusters from scratch, with control over when to read/write from cache
	Inputs:
	V -> |n|x3 geometry
	T -> |T|x4 tet indices
	W -> |n|xm clustering features

	read_cache -> whether or not to attempt to read modes and clusters from cache.
	write_cache -> whether or not to write modes and clusters to cache.
	cache_dir -> directory where clusters cache is
	(Optional)
	recompute_clusters_if_not_found -> whether or not to recompute cluster from scratch if not found in cache (default true)
	*/
	void init_with_cache(MatrixXd& V, MatrixXi& T, MatrixXd& W, bool read_cache, bool write_cache,
		 string& cache_dir,
		bool recompute_if_not_found = true)
	{
		if (read_cache)
		{
			bool found_subspace_cache = this->read_from_cache_recompute(V, T, 
				W, cache_dir, recompute_if_not_found,
				write_cache);
			
		}
		else
		{
			this->init(V, T, W);
			if (write_cache)
				this->write_to_cache(cache_dir);
		}
	}

	/*
	Reads modes, clusters from cache directories. If they are not found, recompute each one of these as specified by
	recomp_modes_if_not_found//
	*/
	bool read_from_cache_recompute(MatrixXd& V, MatrixXi& T, MatrixXd& W,  std::string& clusters_dir,
		bool recomp_if_not_found, bool write_cache)
	{
		bool success = true;

		
		success = success && read_clusters_from_cache(clusters_dir);
		if (!success && recomp_if_not_found)
		{
			printf("Computing %i clusters from scratch... \n",
				num_clusters);
			MatrixXd C;// cluster centers in weight space... just throw this out
				//compute_clusters_displacement_features(T, B, L, params.num_clusters, params.num_clustering_features, l, C, true);
			VectorXd L = VectorXd::Ones(W.cols());

			compute_clusters_weight_features(T, W, L, num_clusters, W.cols(), l, C, false);
		
			if (write_cache)
				write_clusters_to_cache(clusters_dir);
		}
		return success;
	}

	/*
	Writes modes and clusters to cache directories
	clusters_dir -> (string) where to save clusters directory (labels.DMAT)
	*/
	bool write_to_cache( std::string& clusters_dir)
	{
		bool success = true;;
		success = success && write_clusters_to_cache(clusters_dir);
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
	  clusters_dur - (string) directory where cluster labels.DMAT is stored
	*/
	bool read_from_cache( std::string& clusters_dir)
	{
		bool success = true;
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
			bool correct_num_clusters = l.maxCoeff() + 1 < num_clusters + 10 && l.maxCoeff() + 1 > num_clusters - 10;

		success = success && correct_num_clusters;
		
		printf("Successfully Read  Clusters from cache dir %s! \n", clusters_dir.c_str());
	
		return success;
	}

};