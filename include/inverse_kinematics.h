#pragma once
#include "Skeleton.h"

namespace ik
{
enum IK_SOLVER { GD, LM };
/*
Performs inverse kinematic on our skeleton. Finds xzx angle for each bone that enforces boundary conditions of xb0 at indices b.
Implements both gradient descent and Leventberg Marquad (only former operational for now)
Stores the result in place inside input skeleton.
*/
void inverse_kinematics(Skeleton& skeleton, Eigen::VectorXi& b, Eigen::VectorXd& xb0, int max_iter = 100, double tol = 1e-10, IK_SOLVER solver = GD, double mu=1e-5);
}