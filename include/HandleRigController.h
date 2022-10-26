#pragma once
#include <Eigen/Core>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include "RigController.h"
#include "fast_cd_viewer.h"
/*
Class that takes into account general handle rig parameters p0. This is a 12xh vector, containing 12 entries for each handle.
This class links the handle rig to a way of interacting with it using mouth controls. This involes a few things: it draws on the viewer
points for each handle, it lets you select which handle you want to manipulate, and finally it lets you move each handle with an ImGuizmo
*/
class HandleRigMouseController : public RigController
{
public:
	HandleRigMouseController() {};


	HandleRigMouseController(Eigen::VectorXd& p0, fast_cd_viewer* viewer, std::string animation_dir = "");

	HandleRigMouseController(Eigen::VectorXd& p0,  igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo, std::string animation_dir="");
	/*
	returns the constrained displacement of all vertices pinned by this rig.
	*/
	virtual void get_pinned_motion(Eigen::MatrixXd& X, std::vector<Eigen::VectorXi> bI, Eigen::VectorXd& bc);

	virtual void render(igl::opengl::glfw::Viewer& viewer);

	virtual bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return false; };

	virtual bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {return false;};

	virtual bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

	virtual bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier);

	virtual void init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);
	

	Eigen::MatrixX4f matrix4f_from_parameters(Eigen::VectorXd& p, int i);

	virtual void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);

	void reset();


	virtual void set_scripted_motion(int step);
	
	

public:


	Eigen::MatrixXd V; //point handle positions.
	int handleI;
	int vis_id;

	igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo;




};