#pragma once
#include "RigController.h"
#include "fast_cd_init.h"
#include "update_parameters_at_handle.h"
#include <Eigen/Geometry>
/*
Script that moves every handle in rig left and right by the same amount.
*/
class RotatingRigController : public RigController
{

double rad_per_frame;
Eigen::Vector3d rot_axis;
public:
	RotatingRigController(FastCDInit * init)
	{ 
		rad_per_frame = ((InitSimRotation*)init)->rad_per_frame;
		p_rel.resize(12);
	};

	virtual Eigen::VectorXd query_rel(int step)
	{
		double angle = step * rad_per_frame;
		//Get rotation matrix
		Eigen::Matrix3f m = Eigen::AngleAxisd(angle, rot_axis).toRotationMatrix();

		Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
		T.block(0, 3, 3, 1) = m;

		update_parameters_at_handle(p_rel, T.cast<float>(), 0);
	}

	
};

