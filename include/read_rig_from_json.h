#pragma once

#include <fstream>
#include <Eigen/Core>
#include <json.hpp>
#include <filesystem>
#include <igl/list_to_matrix.h>
using namespace nlohmann;
using namespace std;
using namespace Eigen;

/*
Reads rig from json file 
Inputs :
	rig_path - (string) where is my rig.json file?
Outputs:
	W - n x m matrix of rig weights
	P0 - 4mx3 matrix of rest bone world transforms
		, where each 4x3 block is an affine handle
	pI - mx1 parent bone Indicies, in the case of a skeleton rig
	l - mx1 bone lengths, useful for visualization
	V - nx 3 mesh geometry
	F - Fx 3|4 mesh face/tet indices
    rig_type - either "surface" or "volume", 
*/
void read_rig_from_json(string& rig_path, MatrixXd& W, MatrixXd& P0, VectorXi& pI, VectorXd& l, MatrixXd& V, MatrixXi& F, string& rig_type)
{

	namespace fs = std::filesystem;
	if (!fs::exists(fs::path(rig_path)))
	{
		printf("Could not find rig path .json file at %s \n", rig_path.c_str());
		return;
	}

	json j;
	std::ifstream i(rig_path);
	i >> j;

	std::vector<int> pI_list = j.value("pI", std::vector<int>() );
	std::vector<double> l_list = j.value("lengths", std::vector<double>());

	std::vector<std::vector<double>> V_list = j["V"];
	std::vector<std::vector<double>> W_list = j["W"];
	std::vector<std::vector<int>> F_list = j["F"];
	std::vector<std::vector<std::vector<double>>> P0_list = j["p0"];
	int d = V_list[0].size();
	int n = P0_list.size();
	P0.resize(P0_list.size() * (d + 1), 3);

	Matrix3d R;

	for (int i = 0; i < n; i++)
	{
		MatrixXd M;
		igl::list_to_matrix(P0_list[i], M);
		if (M.rows() == d+1)
			P0.block(i * (d + 1), 0, d + 1, d) = M;
		else if (M.rows() == d)
			P0.block(i * (d + 1), 0, d + 1, d) = M.transpose();
	}
	//igl::list_to_matrix(P0_list, P0);
	igl::list_to_matrix(V_list, V);
	igl::list_to_matrix(F_list, F);
	igl::list_to_matrix(W_list, W);
	igl::list_to_matrix(pI_list, pI);
	igl::list_to_matrix(l_list, l);

	rig_type = j.value("rig_type", "surface");
	if (rig_type != "surface" && rig_type != "volume")
	{
		printf("Rig type %s not supported, must be either surface or volume", rig_type.c_str());
	}
	//i.close();
}
void read_rig_from_json(string& rig_path, MatrixXd& W, MatrixXd& P0, MatrixXd& V, MatrixXi& F, string& rig_type)
{
	VectorXi pI;
	VectorXd l;
	read_rig_from_json(rig_path, W, P0, pI, l, V, F, rig_type);

}