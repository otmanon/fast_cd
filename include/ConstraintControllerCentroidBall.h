#pragma once
#include "ConstraintController.h"

class ConstraintControllerCentroidBall : public ConstraintController
{
public:
	ConstraintControllerCentroidBall(Eigen::MatrixXd& X, double r, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);

	void render(igl::opengl::glfw::Viewer& viewer);

	void get_interactive_motion(Eigen::MatrixXd& P);

	void init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo) {};

public:
	double r;
	Eigen::Matrix4f T_rest_inv;
	Eigen::Matrix4f T;
};