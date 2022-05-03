#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

class FastCDViewer
{
public:
	FastCDViewer();

	FastCDViewer(igl::opengl::glfw::Viewer* viewer);

;

	void launch();

	void set_pre_draw_callback(std::function<bool()>& callback);

	void configure_clusters( Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters);

	void configure_deformation_texture( Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W);

	void configure_color_texture(std::string texture_filepath, Eigen::MatrixXd& V_coarse, Eigen::MatrixXi& F_coarse, Eigen::MatrixXd& V_fine, Eigen::MatrixXi& F_fine,
		Eigen::MatrixXd& UV_fine, Eigen::MatrixXi& FUV_fine);

	void configure_matcap(std::string matcap_file, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

	/*
	Renders mesh  full space mesh V.
	*/
	void render_full(Eigen::MatrixXd& V, int id);


	/*
	Renders mesh with reduced deformation z, and rig parameters p.
	Must do a projection step Bz + Jp to get the final V. Does this step on the CPU
	*/
	void render_reduced_cpu_proj(Eigen::VectorXd& z, Eigen::VectorXd & p, Eigen::MatrixXd& B, Eigen::SparseMatrix<double>& J, int id);

	void render_reduced_gpu_proj(Eigen::VectorXd& z, Eigen::VectorXd& p, int n, int id);

	//void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);


public:
	igl::opengl::glfw::Viewer* viewer;

	std::function<bool()> callback;
	int fid, cid; //indices in the viewer.data_list array corresponding to the fine mesh, and the coarse mesh
};