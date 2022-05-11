#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
class FastCDViewer
{
public:
	FastCDViewer();

	FastCDViewer(igl::opengl::glfw::Viewer* viewer);

	//attach guizmo too
	FastCDViewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuiWidget* guizmo);
	void launch();

	void set_pre_draw_callback(std::function<bool()>& callback);

	void attach_guizmo(igl::opengl::glfw::imgui::ImGuizmoWidget& guizmo)
	{
		igl::opengl::glfw::imgui::ImGuiPlugin  * imgui_plugin = new igl::opengl::glfw::imgui::ImGuiPlugin();
		igl_v->plugins.push_back(imgui_plugin);
		// push back menu here
		imgui_plugin->widgets.push_back(&guizmo);
		//imgui_plugin->widgets.push_back(guizmo);
	}

	void set_lighting_factor(double f)
	{
		igl_v->core().lighting_factor = f;
	}
	void clear_all()
	{
		for (int i = 0; i < igl_v->data_list.size(); i++)
		{
			igl_v->data_list[i].clear();
		}
	}
	void set_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int id)
	{
		igl_v->data_list[id].set_mesh(V, F);
	}

	void set_camera_zoom(double zoom)
	{
		igl_v->core().camera_zoom = zoom;
	}
	void set_vertices(Eigen::MatrixXd& V, int id)
	{
		igl_v->data_list[id].set_vertices(V);
	}
	void compute_normals(int id)
	{
		igl_v->data_list[id].compute_normals();
	}

	void clear(int id)
	{
		igl_v->data_list[id].clear();
	}

	void set_color(Eigen::RowVector3d & d, int id)
	{
		igl_v->data_list[id].set_colors(d);
	}
	void set_color(Eigen::MatrixXd& d, int id)
	{
		igl_v->data_list[id].set_colors(d);
	}

	void invert_normals(int id)
	{
		igl_v->data_list[id].invert_normals = !igl_v->data_list[id].invert_normals;
	}

	void set_points(Eigen::MatrixXd& points, int id)
	{
		Eigen::RowVector3d z(0, 0, 0);
		igl_v->data_list[id].set_points(points, z);
	}

	void set_points(Eigen::MatrixXd& points, Eigen::RowVector3d& z, int id)
	{
		igl_v->data_list[id].set_points(points, z);
	}

	void set_points(Eigen::MatrixXd& points, Eigen::MatrixXd& color, int id)
	{
		igl_v->data_list[id].set_points(points, color);
	}

	void set_lines(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXd& color, int id)
	{
		igl_v->data_list[id].set_edges(V, E, color);
	}

	void set_lines(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::RowVector3d& z, int id)
	{
		igl_v->data_list[id].set_edges(V, E, z);
	}

	void set_lines(Eigen::MatrixXd& V, Eigen::MatrixXi& E, int id)
	{
		Eigen::RowVector3d z(0, 0, 0);
		igl_v->data_list[id].set_edges(V, E, z);
	}
	void set_face_based(bool face_based, int id)
	{
		igl_v->data_list[id].face_based = face_based;
	}
	void set_show_lines(bool show_lines, int id)
	{
		igl_v->data_list[id].show_lines = show_lines;
	}
	void set_show_faces(bool show_faces, int id)
	{
		igl_v->data_list[id].show_faces = show_faces;
	}
	void configure_clusters( Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters);

	void configure_clusters(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters, int id);

	void configure_deformation_texture( Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W);

	void configure_solid_color_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::RowVector3d& color, int id);
	void configure_color_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& color, int id); 

	void configure_color_texture(std::string texture_filepath, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
		Eigen::MatrixXd& TC, Eigen::MatrixXi& FTC, int id = 0);

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
	igl::opengl::glfw::Viewer* igl_v;
	igl::opengl::glfw::imgui::ImGuiWidget* guizmo; //add guizmo

	std::function<bool()> callback;
	int fid, cid; //indices in the viewer.data_list array corresponding to the fine mesh, and the coarse mesh
};