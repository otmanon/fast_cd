#include "read_anim_file.h"
#include <iostream>
#include <fstream>
#include <json.hpp>
#include <stdio.h>
#include <stdio.h>
#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
using namespace nlohmann;
std::vector<std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> > read_anim_file( std::string file, const Skeleton& skeleton, Eigen::Matrix4d& S, Eigen::Vector3d& t)
{
	#ifdef WIN32
	namespace fs = std::filesystem;
	#else
	namespace fs = std::experimental::filesystem;
	#endif
	
	std::vector<std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> > anim_T;
	if (!fs::exists(fs::path(file)))
	{
		std::cout << "no animation found" << std::endl;
		return anim_T;
	}
	else
	{
		std::ifstream i(file);
		json j;
		i >> j;

		std::vector<std::vector<std::vector<std::vector<double>>>> list = j["anim_world_space_matrices"];
		anim_T.resize(list.size());
		Eigen::Matrix4d A;
		std::vector<std::vector<double>>a;
		for (int i = 0; i < list.size(); i++)
		{
			anim_T[i].resize(list[i].size());
			for (int bi = 0; bi < list[i].size(); bi++)
			{
				a = list[i][bi];
				A <<	a[0][0], a[0][1], a[0][2], a[0][3],
						a[1][0], a[1][1], a[1][2], a[1][3],
						a[2][0], a[2][1], a[2][2], a[2][3],
						a[3][0], a[3][1], a[3][2], a[3][3];
				anim_T[i][bi].matrix() = A; // blender Z is up, not forward! rotate about x axis
				anim_T[i][bi] = Eigen::AngleAxisd(-igl::PI * 0.5 , Eigen::Vector3d::UnitX()) * anim_T[i][bi];
				//anim_T[i][bi].matrix() = S * anim_T[i][bi].matrix();
				//anim_T[i][bi].matrix().block(0, 3, 3, 1) += t;
				
				//anim_T is a matrix that goes from local bone space to world (in blender)
				// local bone space in blender has bone length along the y axis. our local bone space assumes bone length along the x axis.
				// to accomodate this, we need to prerotate about z by 90 degrees
				anim_T[i][bi] = anim_T[i][bi] * Eigen::AngleAxisd( igl::PI *0.5,  Eigen::Vector3d::UnitZ()) *skeleton[bi].rest_T.inverse();  // get skeleton inverse rest pos
			}
		}

	}
	return anim_T;
}