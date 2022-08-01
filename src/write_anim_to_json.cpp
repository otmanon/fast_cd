#include "write_anim_to_json.h"
#include "flattened_parameters_to_Affine3D_list.h"
#include "format_transforms.h"

#include <json.hpp>
#include <fstream>
#include <filesystem>
bool write_anim_to_json(const Eigen::MatrixXd& P,  const std::string file, std::string const format)
{
	std::vector<std::vector<std::vector<std::vector<double>>>> anim;
	for (int i = 0; i < P.cols(); i++)
	{
		std::vector<Eigen::Affine3d> T, T_format;
		Eigen::VectorXd p = P.col(i);
		
		flattened_parameters_to_Affine3d_list(p, T);
		
		format_transforms(T, format, T_format); //blender rotates about x 90 degrees (flipping z and y)
		
		std::vector<std::vector<std::vector<double>>> frame;
		for (int b = 0; b < T.size(); b++)
		{
			std::vector<std::vector<double>> mat;
			for (int u = 0; u < 3; u++)
			{
				std::vector<double> row;
				for (int v = 0; v < 4; v++)
				{
					row.push_back(T[b](u, v));
				}
				mat.push_back(row);
			}
			frame.push_back(mat);
		}
		anim.push_back(frame);
	}

	namespace fs = std::filesystem;
	if (!fs::exists(fs::path(file).parent_path()))
	{

		fs::create_directories(fs::path(file).parent_path());
	}

	//write rig to json
	std::ofstream o(file);
	nlohmann::json j;
	j["anim_world_space_matrices"] = anim;
	o << std::setw(4) << j << std::endl;

	return true;
}