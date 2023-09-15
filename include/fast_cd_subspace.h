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

/*
This class is in charge of constructing the fast_cd subspace

*/
struct  fast_cd_subspace
{

	int num_modes;
	int num_clusters;
	std::string mode_type;
	MatrixXd B;    // final subspace
	MatrixXd W;    //skinning weights (if applicable)
	VectorXi l;   // clusters

	fast_cd_subspace() {};
	fast_cd_subspace(std::string modes_dir, std::string clusters_dir, std::string mode_type, int num_modes, int num_clusters) {
		this->mode_type = mode_type;
		this->num_clusters = num_clusters;
		this->num_modes = num_modes;
		read_from_cache(modes_dir, clusters_dir);
	};



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

		//allow for some wiggle room in cluster spec
		if (!success)
			printf("labels.DMAT file not found in clusters cache dir specified.");

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

		if (mode_type == "skinning")
		{
			success = success && igl::readDMAT(modes_dir + "/W.DMAT", W);

			if (W.cols() >= num_modes)
			{
				W = W.block(0, 0, W.rows(), num_modes); //
			}
			else
				success = false;

		}
		else if (mode_type == "displacement")
		{
			success = success && igl::readDMAT(modes_dir + "/B.DMAT", B);
			if (B.cols() >= num_modes)
			{
				B = B.block(0, 0, B.rows(), num_modes);
			
			}
			else
				success = false;

		}

		if (!success)
		{
			printf("Coult not find W.DMAT or B.DMAT with a corresponding L.DMAT file in modes cache dir specified");
		}
		return success;
	}
};