#pragma once
#include "HandleRigController.h"

/*
Script that moves every handle in rig left and right by the same amount. 
*/
class LeftRightScriptedHandleRigController : public HandleRigMouseController
{

public:
	LeftRightScriptedHandleRigController(Eigen::VectorXd& p0, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo) :HandleRigMouseController(p0, viewer, guizmo) {};

	void set_scripted_motion(int step);

};