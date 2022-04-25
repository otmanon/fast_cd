#pragma once
#include "Skeleton.h"
/*
Applies affine transformation to each bone in order to get the world space positions. Gets the affine transformation through the xyx euler angles of each bone
*/
void forward_kinematics(const Skeleton& bone_list, std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T);


/*
Applies affine transformation to each bone in order to get the world space positions. Gets the affine transformation through a list of input matrices T_basis

T_basis is a local (to each bone) rotation matrix relative to its rest pose, and it'as parent'as local rotation matrix

T_world = T_parent_world * T_basis * T_local
*/
void forward_kinematics(const Skeleton& bone_list, std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T_basis, std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T);

