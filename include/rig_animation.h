#pragma once

#include <string>
#include <Eigen/Core>

#include <read_anim_file.h>
#include <load_animation.h>
#include <get_relative_parameters.h>
class RigAnimation
{
public:
	RigAnimation(std::string anim_path, Eigen::VectorXd& p0, double scale, Eigen::Vector3d& t)
	{
		bool _n;
	}

	RigAnimation(std::string anim_path, Eigen::VectorXd& p0, Eigen::VectorXi& pI)
	{
		bool _n;
		this->p0 = p0;
		load_animation_and_fit_to_rest_skeleton_json(anim_path, p0, pI, P, _n);
	}

	Eigen::VectorXd query_rel(int step)
	{
		Eigen::VectorXd p_rel;
		Eigen::VectorXd p_glob= P.col(step);
		//  get_relative_parameters(p0,p0, p_rel0);
		get_relative_parameters(p0, p_glob, p_rel);
		return p_rel;
	}

	Eigen::VectorXd query_glob(int step)
	{
	
		Eigen::VectorXd p_glob = P.col(step);
		return p_glob;
	}

	//void fit(std::string rig_path, double scale, Eigen::Vector3d& t)
	//{
	//
	//}
	//
	//void query(int step, Eigen::VectorXd& p_rel)
	//{
	//
	//}
public:

	Eigen::VectorXd p0;
	Eigen::MatrixXd P;
	Eigen::MatrixXd P_rel;

};