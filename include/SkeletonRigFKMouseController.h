#pragma once
#include "HandleRigController.h"


	 /*
	 From absolute world rig parameters p, build your skeleton mesh.
	 */
 void get_skeleton_mesh(float thickness, Eigen::VectorXd& p, Eigen::VectorXd& bl, Eigen::MatrixXd& renderV, Eigen::MatrixXi& renderF, Eigen::MatrixXd& renderC);
class SkeletonRigFKMouseController : public HandleRigMouseController
{

public:
	/*
	Creates a skeleton Rig Mouse Controller, controlling a SkeletonRig using forward kinematics, making use of the skeleton hierarchy

	Inputs:
	p0 = 12bx1 vector of row order flattened rest-pose rig parameters
	pI = bx1 vector of parent indices into our bones (-1 if bone does not have a parent)
	pl = bx1 vector of bone lengths
	viewer = the libigl viewer object
	guizmo = the guizmo, declared outside of this class, but this class re-initializes it

	*/
	 SkeletonRigFKMouseController(Eigen::VectorXd& p0, Eigen::VectorXi& pI , Eigen::VectorXd& pl, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget * guizmo, std::string animation_dir="");


	 //show the current skeleton visually
	 void init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);

	 //show your mesh girl!
	 void render(igl::opengl::glfw::Viewer& viewer);


	 virtual void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);
	 bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
	 void reset();

	 void set_scripted_motion(int step);
public:
	Eigen::VectorXi pI; //parent indices
	Eigen::VectorXd pl; //bone lengths


	//name of the game here is to ultimately end up with rig parameters relative to the rest pose of their respective bones

	//Guizmo gives us the global transformation of each rig parameter... we need to transform that into a transformation relative to bones parent, 
	// then propagate that through our skeleton hierarchy using forward kinematics. This would give us the new global transformation matrix for each bone
	// finally, we need to convert that back fown to a relative transformation, flatten it, and give it to p_rel, which is what is ultimately needed for our rig

	//converting these to matrices first is way easier for our current implementation of FK
	std::vector<Eigen::Affine3d> A0; //list of rest pose rig transofrmation matrices
	std::vector<Eigen::Affine3d> A_rel_parent; // list of transformation matrices relative to parent
	std::vector<Eigen::Affine3d> A_rel; // list of transformation matrices relative to rig

	//override HandleRIg's V to now be bone TIPS, not handle origins
	
	
	std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>R_rel_parent; //relative rotations relative to parent

	Eigen::MatrixXd C;		//positions of all the tips AND the joints
	Eigen::MatrixXd C0;     //ORIGINAL positions of all the tips and the joints (as needed by FK)
	Eigen::MatrixXi BE; // bone edges, indexing joints C

	Eigen::MatrixXd renderV,renderC;
	Eigen::MatrixXi  renderF;

	float thickness;


	bool pause;
	bool loaded_anim;
};