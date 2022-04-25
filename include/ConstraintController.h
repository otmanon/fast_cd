#pragma once
#include <Eigen\Core>
#include <Eigen\Sparse>
#include <igl\opengl\glfw\Viewer.h>
#include <igl\opengl\glfw\imgui\ImGuizmoWidget.h>

class ConstraintController
{
public:
	ConstraintController() {};
	ConstraintController(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo) { 
	//	guizmo->visible = false;
	//	guizmo->callback = [&](const Eigen::Matrix4f& A)
	//	{
//
	//	};
	//	guizmo->T.setIdentity();
		this->viewer = viewer;
		this->guizmo = guizmo;
		apply = false; };

	/*
	Returns the desired boundary condition values of the vertices in question.
	*/
	virtual void get_interactive_motion(Eigen::MatrixXd& P) {};

	/*
	Returns the desired boundary condition values of the vertices in question.
	*/
	virtual void get_scripted_motion(Eigen::MatrixXd& P){};

	/*
	Returns the linear equality constraint matrix
	*/
	virtual void linear_equality_constraint_matrix(Eigen::SparseMatrix<double>& S){};

	/*
	Linear equality constraint rhs
	*/
	virtual void linear_equality_constraint_rhs(Eigen::VectorXd& rhs){};

	/*
	Keycallback that controls
	*/
	virtual bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier) { return false; };

	virtual void reset() { bc = bc0; };

	virtual void init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo) {};
	
	//Give these default functions, even tho we may not have to use all of them. Be sure to call them in the animation interfaces respective callbacks
	virtual void render(igl::opengl::glfw::Viewer& viewer) {};

	virtual bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return false; };

	virtual bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return false; };

	virtual bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return false; };



public:
	bool apply; //are we applying the constraint or not?
	Eigen::MatrixXd bc;
	Eigen::MatrixXd bc0;
	Eigen::VectorXi bI;

	Eigen::SparseMatrix<double> S;
	igl::opengl::glfw::Viewer* viewer; 
	igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo;
};