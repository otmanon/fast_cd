#pragma once
#include <Eigen/Core>
#include <fstream>
#include <json.hpp>
#include <igl/list_to_matrix.h>
#include <filesystem>
using namespace nlohmann;
using namespace std;
using namespace Eigen;

/*
Reads a rig animation from a json file.
Input:
	anim_path - path to anim.json file
Output:
	animP - 12m x #frames matrix containing the world 
				space animations bone transformations 
				for each frame 
*/
void read_rig_anim_from_json(string& anim_path, MatrixXd& animP)
{
	json j;
	std::ifstream i(anim_path);
	i >> j;

	namespace fs = std::filesystem;
	if (!fs::exists(fs::path(anim_path)))
	{
		printf("Could not find animation .json file at %s \n", anim_path.c_str());
		return;
	}
	std::vector<std::vector<std::vector<std::vector<double>>>> animP_list = j["P"];
	
	assert(animP_list[0][0].size() == 3 && animP_list[0][0][0].size() == 4);
	int nf = animP_list.size(); // num frames
	int nb = animP_list[0].size();

	MatrixXd Pf;
	Pf.resize(4 * nb, 3);
	VectorXd pf;
	MatrixXd M;
	animP.resize(nb * 12, nf);
	animP.setZero();
	for (int i = 0; i < nf; i++)
	{
		Pf.setZero();
		for (int b = 0; b < nb; b++)
		{
			igl::list_to_matrix(animP_list[i][b], M);
			Pf.block(b *4, 0, 4, 3) = M.transpose();
		}
		pf = Map<VectorXd>(Pf.data(), 12 * nb, 1);
		animP.col(i) = pf;
	}
	
}