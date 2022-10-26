#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
class fast_cd_viewer
{
public:
	fast_cd_viewer();

	fast_cd_viewer(igl::opengl::glfw::Viewer* viewer);

	//attach guizmo too
	fast_cd_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);
	void launch();


	bool default_key_pressed_callback(igl::opengl::glfw::Viewer& viewer, unsigned int unicode_key, int modifiers, int id);

	void set_pre_draw_callback(std::function<void()>& callback);

	void set_key_pressed_callback(std::function<bool(unsigned int, int)>& callback_key_pressed);

	void set_mouse_down_callback(std::function<void(int button, int modifier)>& callback_mouse_down) {
		igl_v->callback_mouse_down = [&](igl::opengl::glfw::Viewer& v,  int button, int modifier)->bool
		{
			callback_mouse_down(button, modifier);
			return false;
		};
	};

	void set_mouse_down_callback(std::function<bool(igl::opengl::glfw::Viewer& viewer, int button, int modifier)>& callback_mouse_down) {
		igl_v->callback_mouse_down = callback_mouse_down;
	};

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

	void set_camera_eye(Eigen::RowVector3d& p)
	{
		igl_v->core().camera_eye = p.cast<float>();
	}

	void set_camera_center(Eigen::RowVector3d& c)
	{
		igl_v->core().camera_center = c.cast<float>();
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

	void set_double_sided(bool ds, int id)
	{
		igl_v->data_list[id].double_sided = ds;
	}
	void invert_normals(bool invert_normals, int id)
	{
		igl_v->data_list[id].invert_normals = invert_normals;
	}

	void set_points(Eigen::MatrixXd& points, int id)
	{
		Eigen::RowVector3d z(0, 0, 0);
		igl_v->data_list[id].set_points(points, z);
	}

	void configure_wireframe_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int id)
	{
		igl_v->data_list[id].clear();
		igl_v->data_list[id].set_mesh(V, F);
		igl_v->data_list[id].show_faces = false;
		igl_v->data_list[id].show_lines = true;
	}

	void set_data(Eigen::VectorXd &D, int i)
	{
		igl_v->data_list[i].set_data(D);
	}

	void set_data(Eigen::VectorXd& D,  igl::ColorMapType cmap, int i)
	{
		igl_v->data_list[i].set_data(D, cmap);
	}

	void set_data(Eigen::VectorXd& D, double min, double max, int num_steps, igl::ColorMapType cmap, int i)
	{
		igl_v->data_list[i].set_data(D, min, max, cmap, num_steps );
	}

	void set_data(Eigen::VectorXd& D, double min, double max, int i)
	{
		igl_v->data_list[i].set_data(D, min, max);
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


	void set_data_colormap(Eigen::VectorXd& d, Eigen::MatrixXd& cmap, int id)
	{
		igl_v->data_list[id].set_data(d);
		igl_v->data_list[id].set_colormap(cmap);
	}





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

	/*
	GUIZMO controls
	*/
	void init_guizmo(bool visible, Eigen::Matrix4f A0, std::function<void(const Eigen::Matrix4f& A)>& callback, ImGuizmo::OPERATION = ImGuizmo::TRANSLATE);

	void set_guizmo_operation(ImGuizmo::OPERATION op)
	{
		guizmo->operation = op;
	};
	ImGuizmo::OPERATION get_guizmo_operation()
	{
		return guizmo->operation;
	}
public:
	igl::opengl::glfw::Viewer* igl_v;
	igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo; //add guizmo
	igl::opengl::glfw::imgui::ImGuiPlugin* imgui_plugin;

	igl::opengl::glfw::imgui::ImGuiMenu* imgui_menu;
	std::function<bool()> callback;
	int fid, cid; //indices in the viewer.data_list array corresponding to the fine mesh, and the coarse mesh
};