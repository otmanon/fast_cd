#pragma once
#include <Eigen\Core>
/*
Creates a rig for a mesh that has two handles, one on the lft most and one on the right most part of the mesh.

These two handles only control the left most/right most part of the mesh. All those weights are 1. All interior weights are 0.
*/

void create_two_handle_rig(std::string file_path, Eigen::MatrixXd& V);