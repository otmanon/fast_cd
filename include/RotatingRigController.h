#pragma once
#include "RigController.h"
#include "fast_cd_init.h"
#include "update_parameters_at_handle.h"
#include <Eigen/Geometry>
#include "ScriptedRigController.h"
/*
Script that assumes an affine rig and rotates the first bone 90 degrees
*/
class Rotater : public ScriptedRigController
{
double rad_per_frame;
Eigen::Vector3d rot_axis;
public:
	Rotater(Eigen::VectorXd p0, FastCDInit * init)
	{ 
		rad_per_frame = ((InitSimScriptedRig*)init)->rad_per_frame;
		p_rel.resize(p0.rows());
		max_step = ((InitSimScriptedRig*)init)->max_step;
		rot_axis = ((InitSimScriptedRig*)init)->rot_axis;
	};

	virtual Eigen::VectorXd query_rel(int step)
	{
		double angle = step * rad_per_frame;
		//Get rotation matrix
		Eigen::Matrix3d m = Eigen::AngleAxisd(angle, rot_axis).toRotationMatrix();

		Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
		T.block(0, 0, 3, 3) = m;

		
		update_parameters_at_handle(p_rel, T.cast<float>(), 0);
		return p_rel;
	}

	
};

