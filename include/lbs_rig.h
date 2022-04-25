#pragma once
#include "rig.h"
#include <json.hpp>
#include "inverse_kinematics.h"
#include <igl/opengl/glfw/Viewer.h>
using json = nlohmann::json;

enum INTERACTION_MODE{INVERSE_KINEMATICS, ANIMATION };
class LBSRig : public Rig
{
public:
	LBSRig() {};

	LBSRig(std::string filename, Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::VectorXi& bI, std::string anim_file_dir = "");

	/*
	Given our surface weights, performs #B laplace solves to find correspondibg set of #B volume weights
	Input:
	surfaceW - |surfaceX|x|B| matrix of weights, how each bone column affect each surface vertex row
	bI -			list of indices in X corresponding to surface indices
	X -				list of full vertex positions
	T -				tet mesh faces indexing X
	*/
	Eigen::MatrixXd surface_to_volume_weights(Eigen::MatrixXd& surfaceW, Eigen::VectorXi& bI, Eigen::MatrixXd& X, Eigen::MatrixXi& T);

	Eigen::MatrixXd read_vertices_from_json(json& j);

	Eigen::MatrixXi read_faces_from_json(json& j);

	Skeleton read_bones_from_json(json& j);

	Eigen::MatrixXd weight_matrix_from_bones(Skeleton& bone);

	Eigen::MatrixXd get_all_bone_heads(Skeleton& bone);


	//will soon be removed... cant wait till i can remove this one
	void attach_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin) {};

	void update_rig_parameters(std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T);

	//these functions below here should all be made virtual

	void init_parameters();

	void init_jacobian();

	void get_rig_motion(Eigen::VectorXd& ur);

	void reset();

	void render(igl::opengl::glfw::Viewer& viewer);

	void init_viewer(igl::opengl::glfw::Viewer& viewer);


	void poll_rig_changes();

	void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);

	bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

	bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

	bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
	
	/*
	Sometimes we renormalize X so that it'as unit length or whatever scale, and translate it so that it'as centered at the origin
	The bones however are of length that scale with surfaceX
	This method figures out the scaling transform S such that  S surfaceX' = X' and applies it to the rest_T of each bone.
	*/
	void fit_bone_transforms_to_input_X();

	void fit_anim_transforms_to_input_X();
public:

	std::vector<std::string> anim_paths, anim_names;
	int current_anim_id, new_anim_id;
	//TODO: get rid of surfaceF asap
	Eigen::MatrixXi surfaceF;	//mesh faces... indexing surfaceX . Only need this so long as we use this surfaceF to draw in the simulationHook
	Eigen::MatrixXd surfaceX;	//mesh surface X
	Eigen::VectorXi bI;			//indices into X of surface X. In numpy: surfaceX = X[bI, :]
	Eigen::MatrixXd A;
	Eigen::MatrixXf rest_T; //starting affine transformation matrix. 
	Eigen::MatrixXf rest_T_inv; //starting affine transformation matrix. (inverted)


	Eigen::Matrix4d best_fit_scale;
	Eigen::Vector3d best_fit_translation;

	Skeleton skeleton;

	//weight matrix |srfaceV|x|B and |volumeV|x|B  indicating the effect of each bone on each vertex position
	Eigen::MatrixXd surfaceW, W;

	Eigen::VectorXd xb0;		//flattened handle positions used in simulation
	Eigen::VectorXi b;			//bone indices that are constrained used in simulation

	int num_b;					//number of bones
	int num_p;					//number of rig parameters (num_b x 12), for an affine bone skeleton

	//heads of each bone, used to set control points etc
	Eigen::MatrixXd P;

	//reduced_solver parameters
	INTERACTION_MODE interaction_mode;
	std::vector<std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>> anim_T; // transformation matrices going from rest pose to world pose for each bone.(not relative to parent anymore)
	INTERACTION_MODE new_interaction_mode;
	int frame_i;
	int total_frames;

	ik::IK_SOLVER ik_solver;
	int max_iter;
	int max_line_search_steps;
	float tol;
	float mu;
	int tol_exp;
	////////////////////////VIEWER Related QUANTITIES ////////////////////
	//rendering quantities
	Eigen::MatrixXd renderV;
	Eigen::MatrixXi renderF;
	Eigen::MatrixXd renderC;

	
	//how thick to make the skeleton.
	double thickness;

	Eigen::Vector3d mouse_win, mouse_world;
	Eigen::Vector3d mouse_drag_win, mouse_drag_world;
	int currentHandle;          //which handle index into V are we currently moving
	bool mouse_dragging;        //is mouse dragging
	bool added_handle;			//did we add a handle this step
	bool moved_handle;

	//gizmo attachment
	igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin;

	Eigen::VectorXi handleI;		//which bone index is constrained? update in render but done update xb0 until sim
	Eigen::MatrixXd handleV;
	Eigen::MatrixXd handleC;		
	
	//handle colors. blue normally, turns red if selected.
	Eigen::RowVector3d red = Eigen::RowVector3d(0.8, 0.2, 0.2);
	Eigen::RowVector3d blue = Eigen::RowVector3d(0.2, 0.2, 0.6);

};