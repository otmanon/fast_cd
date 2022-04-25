#pragma once
#include "ConstraintController.h"
class ConstraintControllerSquashStretch : public ConstraintController
{
public:
	ConstraintControllerSquashStretch(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo) : ConstraintController(viewer, guizmo) 
	{};

	ConstraintControllerSquashStretch(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo, float scale = 1e-2);
	
	/*
	Returns the desired boundary condition values of the vertices in question.
	*/
	void get_interactive_motion(Eigen::MatrixXd& P);

	/*
	Keycallback that controls the squash/stretch. Press A/a to squash and D/d to stretch
	*/
	bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier);

public:

	Eigen::Vector3d centroid;

	float scale;
	float squash_factor;

};