#pragma once
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
/*
Parent Class that provides a libigl-viewer interface to interact with rigs. There are many ways of interacting with rigs, including keyboard controllers,
mouse controllers, video-stream controllers, scripted controllers, recorded animation controllers. Each one may specialize to a specific type of rig, so the user must make sure not to employ 
a specialized controller on a rig it is not suited to.
*/
class RigController
{
public:
	RigController() {};

	virtual void render(igl::opengl::glfw::Viewer& viewer) {};

	virtual bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {return false;};

	virtual bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier){return false;};

	virtual bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier){return false;};
	
	virtual bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier){return false;};

	virtual void init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo){};


	virtual void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu) {};
	/*
	returns the constrained displacement of all vertices pinned by this rig.
	*/
	virtual void get_pinned_motion(Eigen::MatrixXd& X, std::vector<Eigen::VectorXi> bI, Eigen::VectorXd& bc) {};

	
	virtual void reset() {};

	/*
	Gets the scripted constraints of indices bI , with rest pose X.
	*/
	virtual void set_scripted_motion(int step) {};

public:
	// rig parameters.row order flatteneing of a 3x4 affine matrices for each handle / joint
	//These rig parameters are row order flattened, as most of the building of the jacobian matrix depends on this assumption.
	Eigen::VectorXd p_rest; //  world-to-rig matrices, mapping a vector in the world space to the joint's rest frame
	Eigen::VectorXd p_rel;  //  rig parameters relative to rest frame of joint. To get full params do P_full =  P_rel * P_rest;
	Eigen::VectorXd p_rest_inv;  //inverse rest frame, useful to figure out P_rel from P_full. P_rel = P_full * P_rest_inv
};