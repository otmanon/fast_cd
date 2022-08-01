#pragma once
#include "rig.h"
#include <json.hpp>
#include "inverse_kinematics.h"
#include <igl/opengl/glfw/Viewer.h>
	using json = nlohmann::json;

//This rig class should not itself be doing updating of rig parameters.
class HandleRig : public Rig
{
	public:
		HandleRig() {};


		HandleRig(std::string filename, double radius = 5e-2);

		HandleRig(std::string filename, Eigen::MatrixXd& W, double radius = 5e-2);
		/*
		Build single handle at centroid of mesh
		*/
		HandleRig(Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius = 5e-2);

		HandleRig(Eigen::MatrixXd& X, Eigen::VectorXd& p, Eigen::MatrixXd& W, double radius = 5e-2);

		HandleRig(Eigen::MatrixXd& X, std::vector<Eigen::Matrix4f>& P, Eigen::MatrixXd& W, double radius = 5e-2);

		virtual void init_rig_jacobian();

		virtual void get_rig_jacobian(Eigen::SparseMatrix<double>& J);

		virtual void init_rig_selection_matrix(double radius = 0.05);
		void init_null_space();
		//TODO: distinguish between get_rig_params. Get rig motion goes from the input params and should output a full dimensional motion
		virtual void get_rig_motion(Eigen::VectorXd& p, Eigen::VectorXd& ur);


		virtual bool write_rig_to_json(std::string filename);		//stores rest state rig parameters, as well as geometry (X and T), and weight indices (W) to a JSON file. 
		
		virtual bool read_rig_from_json(std::string filename);

		virtual bool write_surface_rig_to_json(std::string filename, Eigen::MatrixXi& T);
	public:

		std::string rig_pinning;

	
};