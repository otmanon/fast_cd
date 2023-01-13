#pragma once

#include <fstream>
#include <Eigen/Core>
#include <json.hpp>
#include <fstream>
#include <filesystem>
#include <igl/matrix_to_list.h>
using namespace nlohmann;
using namespace std;
using namespace Eigen;


bool write_rig_to_json(string& rig_path, const  MatrixXd& W, const MatrixXd& P0, 
	VectorXi& pI, VectorXd& l, const MatrixXd& V, const MatrixXi& T, string& rig_type)
{


	//write rig to json
	std::ofstream o(rig_path);
	nlohmann::json j;

	int d = V.cols();
	std::vector<std::vector<std::vector<double>>> P0_list;
	for (int i = 0; i < P0.rows()/4; i++)
	{
		MatrixXd Pi = P0.block((d+1)* i, 0, (d+1), d);
		std::vector<std::vector<double>> Pi_list;
		igl::matrix_to_list(Pi, Pi_list);
		P0_list.push_back(Pi_list);
	}
	std::vector<std::vector<double>> W_list, V_list;
	std::vector<std::vector<int>> T_list;
	std::vector<double> l_list;
	std::vector<int> pI_list;
	
	igl::matrix_to_list(l, l_list);
	igl::matrix_to_list(pI, pI_list);
	igl::matrix_to_list(T, T_list);
	igl::matrix_to_list(V, V_list);
	igl::matrix_to_list(W, W_list);
	j["rig_type"] = rig_type;
	j["V"] = V_list;
	j["F"] = T_list;
	j["W"] = W_list;
	j["p0"] = P0_list;
	j["lengths"] = l_list;
	j["pI"] = pI_list;

	o << std::setw(4) << j << std::endl;
	o.close();
	return true;
}