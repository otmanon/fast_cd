#pragma once
#include <Eigen/Geometry>
#include <vector>

/*
Assumes p is a stacked list of row order flattened vertices. Assumes 12 parameters per bone
p = [ p1  .......  p2 ........ p3]
Where p1, p2 , p3 are the 12 rig parameters associated with every bone
*/
void flattened_parameters_to_Affine3d_list(Eigen::VectorXd& p, std::vector<Eigen::Affine3d>& A);
