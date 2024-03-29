#pragma once
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <get_absolute_parameters.h>
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

	virtual Eigen::VectorXd query_rel()
	{
		return this->p_rel;
	}

	virtual Eigen::VectorXd query_rel(int step) { return Eigen::VectorXd::Zero(0); };

	virtual Eigen::VectorXd query_rel(double ss, Eigen::RowVector3d& t)
	{
		return this->p_rel;
	}
	virtual Eigen::VectorXd query_glob()
	{
		Eigen::VectorXd p_glob;
		get_absolute_parameters(this->p_rest, this->p_rel, p_glob);
		return p_glob;
	}
public:
	// rig parameters.row order flatteneing of a 3x4 affine matrices for each handle / joint
	//These rig parameters are row order flattened, as most of the building of the jacobian matrix depends on this assumption.
	Eigen::VectorXd p_rest; //  world-to-rig matrices, mapping a vector in the world space to the joint's rest frame
	Eigen::VectorXd p_rel;  //  rig parameters relative to rest frame of joint. To get full params do P_full =  P_rel * P_rest;
	Eigen::VectorXd p_rest_inv;  //inverse rest frame, useful to figure out P_rel from P_full. P_rel = P_full * P_rest_inv


	int current_animation_id;
	std::string animation_dir;
	std::vector<std::string> animation_filepaths;
	std::vector<std::string> animation_filenames;
	Eigen::MatrixXd anim_P;
	bool is_global_anim;			//if true, then anim_P are GLOBAL transformation matrices, not relative as CD is used to. Each parameter needs to then be converted to a relative anim. Mixamo transformations are global


	bool pause;
	bool loaded_anim;
	int anim_step;

	char custom_anim_name[128] = "custom_anim";
	bool recording;

	Eigen::MatrixXd record_P;
};