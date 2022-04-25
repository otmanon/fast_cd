#pragma once
#include "Skeleton.h"

/*
Computes visualizaiton mesh, vertices, faces and colors!
*/
void skeleton_visualization_mesh(
    const Skeleton bone_list,
    const double thickness,
    Eigen::MatrixXd& SV,
    Eigen::MatrixXi& SF,
    Eigen::MatrixXd& SC);

void skeleton_visualization_mesh(
    const Skeleton bone_list,
    std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& frame_T,
    const double thickness,
    Eigen::MatrixXd& SV,
    Eigen::MatrixXi& SF,
    Eigen::MatrixXd& SC);