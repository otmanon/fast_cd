#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>

enum RIG_TYPE { NULL_RIG, AFFINE_RIG, LBS_RIG, HANDLE_RIG };
class Rig
{
	public:
		Rig() {  };

		//Standard setter
		virtual void set_p(Eigen::VectorXd& p) { this->p = p; };

		//initializes rig parameters
		virtual void init_parameters() {};
		//initializes rig jacobian
		virtual void init_jacobian(){};
		//query current rig displacement for each vertex
		virtual void get_rig_motion(Eigen::VectorXd& ur){};
		//optionally intiizializes a gizmo, that can control a simple affine transformation interactively with a mouse 
		virtual void init_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin) {
			this->gizmoPlugin = gizmoPlugin;
			gizmoPlugin->visible = false;
			gizmoPlugin->T = Eigen::Matrix4f::Identity();
			gizmoPlugin->callback = [&](const Eigen::Matrix4f& T)
			{
			};
		};
		
		//be sure to kill gizmo if you delete this object. TODO: Ask libigl admin why I can't kill gizmo mid sim
		virtual void killGizmo(igl::opengl::glfw::imgui::ImGuiPlugin& imgui_plugin) { 
			//delete at front where gizmo should be
			imgui_plugin.widgets.erase(imgui_plugin.widgets.begin());
		};

		//resets rig to its default state
		virtual void reset() {};
		//any rig specific rendering that needs to be done (ie a skeleton, point handles?)
		virtual void render(igl::opengl::glfw::Viewer& viewer){};
		//initializes the viewer for the rig
		virtual void init_viewer(igl::opengl::glfw::Viewer& viewer)
		{ 
			if (viewer.data_list.size() >= vis_id + 1)
				viewer.data_list[vis_id].clear();
		};
		//checks in simulation thread if anything has caused the rig to change in the render thread. If it has, updates simulation accordingly
		virtual void poll_rig_changes(){};
		
		virtual void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu) { };
		//mouse up callback
		virtual bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return false; };
		//mouse move callback
		virtual bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier){return false;};
		//mouse down callback
		virtual bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier){return false;};
		//key callback
		virtual bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier) {return false;};

		//getters for rig jacobian and rig parameters... not necessary if children inherit publically
		Eigen::SparseMatrix<double> get_J()
		{
			return J;
		}

		Eigen::VectorXd get_p()
		{
			return p;
		}


public:
	// rest positions
	Eigen::MatrixXd X;

	//dim
	int dim;

	//number of nodes
	int n;

	//rig parameters
	Eigen::VectorXd p;

	//rest state rig params
	Eigen::VectorXd p0;

	//Weight matrix
	Eigen::MatrixXd W;

	//dim 3|V| x 12sol_p rig jacobian matrix dx/dp... keep as dense for now, mighthave to change it to sparse later
	Eigen::SparseMatrix<double> J;

	//dim 3|V|x3|V| matrix representing the null map of the rig, that is , 1 on the diagonal if the vertex in question is NOT controlled by the rig
	Eigen::SparseMatrix<double> N; 

	//Matrix that selects out vertices bounded by the rig, assuming a radius of mesh width 0.05*mesh_width
	float radius;

	//S u -> flattened displacement of only the constrained vertices
	Eigen::SparseMatrix<double> S;

	//list of vertices attached to each bone. Use this to specify pinned vertex motion.
	std::vector<Eigen::VectorXi> bI;
	//pointer to the viewer because it'as way more convenient if everyone is initialized with this
	igl::opengl::glfw::Viewer* viewer;

	//gizmo attachment. Don't have this for now until I can figure out how to remove a gizmo from widget_list
	igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin;

	int vis_id = 1;

	std::string rig_type;

};

class AffineRig : public Rig
{
public:
	AffineRig() {};

	
	AffineRig(Eigen::MatrixXd& x);

	void init_parameters();

	void init_jacobian();
	
	void get_rig_motion(Eigen::VectorXd& ur);

	void init_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin);




	//Gonna get rid of attach_gizmo soon ok
	void attach_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin);

	void reset();

	bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier);

	void init_null_space();

public:


	Eigen::MatrixXd A;
	Eigen::MatrixXf rest_T; //starting affine transformation matrix. 
	Eigen::MatrixXf rest_T_inv; //starting affine transformation matrix. (inverted)

};