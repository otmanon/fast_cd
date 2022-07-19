#pragma once

#include "HandleRig.h"

void read_bones_from_json(json& j, int num_X, Eigen::MatrixXd& surfaceW, Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::VectorXd& lengths);

/*
opens the rig.json file path provided, reads it and returns 
"surface" -> If the rig only has surface data (as you would get from tools like blender)
"volume" -> If the rig has volumetric data. This either checks the json parameter "format", but if that doesn't exist, it returns volume if the "faces" section has 3 colums.
"null" -> If the rig file does not specify a format, and one cannot be obtained, we return null
*/
std::string get_rig_file_format(std::string filename);

	

class SkeletonRig : public HandleRig
{

public:

//	SkeletonRig(std::string rig, std::string weights, std::string weights_type, Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius = 5e-2);
	/*
	Constructor for an input source file with the surface_file format. This specific case is common and needs to be configured/converted
	to a Volume rig, by diffusing the surface weights. into the mesh.
	*/
	SkeletonRig(std::string surface_file_name, Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius=5e-2);

	SkeletonRig(std::string surface_file_name, Eigen::MatrixXd& W, Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius = 5e-2);

	SkeletonRig(Eigen::MatrixXd& X, Eigen::VectorXd& p0, Eigen::MatrixXd& W, Eigen::VectorXi& pI, Eigen::VectorXd& lengths, double radius = 5e-2);

		//new SkeletonRig(V0, rig->p0, permuted_W, rig->pI, rig->lengths);)

		/*
		Constructor with an input rig file in the volume format. This format just includes mesh vertex rest positions, rest pose rig parameters p0,
		A weight for each bone, bone lengths, and a bone hierarchy list, which is an index into the bone list for each parent
	*/
	SkeletonRig(std::string volume_file_name, double radius=5e-2);

	void init_rig_selection_matrix(double radius = 0.05);

	bool write_rig_to_json(std::string filename);		


	bool read_rig_from_json(std::string filename);


public:
	//length of each bone, used to compute tip
	Eigen::VectorXd lengths;

	//parent index of each bone. 
	Eigen::VectorXi pI;
};