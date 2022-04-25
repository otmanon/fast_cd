#pragma once
#include <Eigen/Core>
#include <igl\unproject_in_mesh.h>

void unproject_in_mesh(const Eigen::Vector2f& win, const Eigen::Matrix4f& view, const Eigen::Matrix4f& proj, const Eigen::Vector4f& viewport, 
	const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::Vector3d& p, std::vector<igl::Hit>& hits);
