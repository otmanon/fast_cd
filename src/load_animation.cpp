#include "load_animation.h"
#include "matrix4f_from_parameters.h"
#include "compute_handle_positions_from_rig_parameters.h"
#include "get_tip_positions_from_parameters.h"
#include "get_joint_positions_from_parameters.h"

#include <iostream>
#include <fstream>
#include <json.hpp>
#include <stdio.h>
#include <stdio.h>

#include "get_bone_length_scales.h"
#include <filesystem>
#include <Eigen/Geometry>
#include <igl/PI.h>
#include <igl/cat.h>
#include <igl/colon.h>
#include <igl/edge_lengths.h>
using namespace nlohmann;

void load_animation(std::string anim_filepath, Eigen::MatrixXd& anim_P, bool is_global)
{

}


void load_animation_and_fit(std::string anim_filepath, Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::MatrixXd& anim_P, bool is_global)
{
	namespace fs = std::filesystem;

	if (!fs::exists(fs::path(anim_filepath)))
	{
		std::cout << "no animation found" << std::endl;
		return;
	}
	else
	{
		std::ifstream i(anim_filepath);
		json j;
		i >> j;

		std::string bone_aligned_to = j.value("bone_aligned_to:", "z");   //blender bones go locally forward in the x direciton i believe
		std::string world_up = j.value("up:", "z"); //blender uses a different world up than us by default 
		std::vector<std::vector<std::vector<std::vector<double>>>> list = j["anim_world_space_matrices"];
		
		anim_P.resize(list[0][0].size()*list[0].size(), list.size());
		Eigen::Matrix4d A, A_ref;
		Eigen::Affine3d T, T_ref;
		std::vector<std::vector<double>> a ;


		
		for (int i = 0; i < list.size(); i++)
		{
			for (int bi = 0; bi < list[i].size(); bi++)
			{
				a = list[i][bi];
				A << a[0][0], a[0][1], a[0][2], a[0][3],
					a[1][0], a[1][1], a[1][2], a[1][3],
					a[2][0], a[2][1], a[2][2], a[2][3],
					a[3][0], a[3][1], a[3][2], a[3][3];
				T.matrix() = A;


					//z is up... change that to y TODO: make thes econditional
				T = Eigen::AngleAxisd(-igl::PI * 0.5, Eigen::Vector3d::UnitX())* T;
		
				// bone is on y, we need it on x
				T = T * Eigen::AngleAxisd(igl::PI * 0.5, Eigen::Vector3d::UnitZ());
				
				//go to bone entry in parameter list
				anim_P.block(12 * bi, i, 4, 1) = T.matrix().row(0);
				anim_P.block(12 * bi + 4, i, 4, 1) = T.matrix().row(1);
				anim_P.block(12 * bi + 8, i, 4, 1) = T.matrix().row(2);

			}
		}

		int num_b = p0.rows() / 12;
		//edge connectivity
		Eigen::VectorXi I1, I2;
		igl::colon(0, num_b - 1, I1);
		I2 = I1.array() + num_b;
		Eigen::MatrixXi BE(num_b, 2);
		BE.col(0) = I1;
		BE.col(1) = I2;

		Eigen::VectorXd p_tmp = anim_P.col(0);
		//anim_P now is all filled up... but entries need to be scaled down
		Eigen::MatrixXd joints;

	
		Eigen::VectorXd lengths;
		//get bone lengths
		double l = get_bone_length_scales(p_tmp, pI);

		//get_tip_positions_from_parameters(p_tmp, pl, tips); //damn Im done... we actually can't do this because pl is just the bone length... need to loop through each bone, go to its parent, calc the diff.
		double l_ref = get_bone_length_scales(p0, pI);

		double s = l / l_ref;

		anim_P /= s;
		

	}
}