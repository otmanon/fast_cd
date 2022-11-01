#pragma once

#include <fstream>
#include <Eigen/Core>
#include <json.hpp>

#include <igl/list_to_matrix.h>
using namespace nlohmann;
using namespace std;
using namespace Eigen;
void read_rig_from_json(string& rig_path, MatrixXd& W, MatrixXd& P0, MatrixXd& V, MatrixXi& F, string& rig_type)
{
	json j;
	std::ifstream i(rig_path);
	i >> j;


	std::vector<std::vector<double>> V_list = j["V"];
	std::vector<std::vector<double>> W_list = j["W"];
	std::vector<std::vector<int>> F_list = j["F"];
	std::vector<std::vector<std::vector<double>>> P0_list = j["p0"];
	int d = P0_list[0].size();
	int n = P0_list.size();
	P0.resize(P0_list.size() * (d+1), 3);

	Matrix3d R;
	
	for (int i = 0; i < n; i++)
	{
		MatrixXd M, M2;
		igl::list_to_matrix(P0_list[i], M);
		
		P0.block(i * (d + 1), 0, d+1, d) = M.transpose();
	}
	//igl::list_to_matrix(P0_list, P0);
	igl::list_to_matrix(V_list, V);
	igl::list_to_matrix(F_list, F);
	igl::list_to_matrix(W_list, W);

	rig_type = j.value("rig_type", "handle");
	//p0 = Eigen::Map<Eigen::VectorXd>(&p0_list[0], p0_list.size());

}