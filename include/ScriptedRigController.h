#pragma once
#include "RigController.h"

class ScriptedRigController : public RigController
{
public:
	ScriptedRigController() {};
	int max_step;
	virtual Eigen::VectorXd query_rel(int step)
	{
		return Eigen::VectorXd::Zero(0);
	};
};