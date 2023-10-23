#pragma once
#include <fast_cd_viewer_parameters.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/png/readPNG.h>

using namespace std;
using namespace Eigen;
class fast_cd_viewer
{
public:
	fast_cd_viewer();

	virtual void  set_menu_callback(std::function<void()>& menu_callback)
	{
		imgui_menu->callback_draw_viewer_menu = [&]()
		{
			menu_callback();
		};
	};

	fast_cd_viewer(igl::opengl::glfw::Viewer* viewer);

	//attach guizmo too
	fast_cd_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);

	virtual void  launch();


	virtual bool default_key_pressed_callback(igl::opengl::glfw::Viewer& viewer, unsigned int unicode_key, int modifiers, int id);

	virtual void set_pre_draw_callback(std::function<void()>& callback);

	virtual void set_key_pressed_callback(std::function<bool(unsigned int, int)>& callback_key_pressed)
	{
        igl_v->callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int modifier)->bool
        {
            printf("Key Pressed");
            return callback_key_pressed(key, modifier);
        };
    }

	virtual void set_mouse_down_callback(std::function<void(int button, int modifier)>& callback_mouse_down) {
		igl_v->callback_mouse_down = [&](igl::opengl::glfw::Viewer& v,  int button, int modifier)->bool
		{
			callback_mouse_down(button, modifier);
			return false;
		};
	};

	virtual void set_mouse_up_callback(std::function<void(int button, int modifier)>& callback_mouse_up) {
		igl_v->callback_mouse_up = [&](igl::opengl::glfw::Viewer& v, int button, int modifier)->bool
		{
			callback_mouse_up(button, modifier);
			return false;
		};
	};

	virtual void set_mouse_down_callback(std::function<bool(igl::opengl::glfw::Viewer& viewer, int button, int modifier)>& callback_mouse_down) {
		igl_v->callback_mouse_down = callback_mouse_down;
	};

	virtual void attach_guizmo(igl::opengl::glfw::imgui::ImGuizmoWidget& guizmo)
	{
		igl::opengl::glfw::imgui::ImGuiPlugin  * imgui_plugin = new igl::opengl::glfw::imgui::ImGuiPlugin();
		igl_v->plugins.push_back(imgui_plugin);
		// push back menu here
		imgui_plugin->widgets.push_back(&guizmo);
		//imgui_plugin->widgets.push_back(guizmo);
	}

	virtual void set_lighting_factor(double f)
	{
		igl_v->core().lighting_factor = f;
	}

	virtual void add_mesh( int& id)
	{
		igl_v->append_mesh();
		id = igl_v->data_list.size() - 1;
		igl_v->data_list[id].clear();
	}
	virtual void add_mesh(const MatrixXd& V, const MatrixXi& F, int& id)
	{
		igl_v->append_mesh();
		id = igl_v->data_list.size() - 1;
		igl_v->data_list[id].clear();
		igl_v->data_list[id].set_mesh(V, F);
	}

	virtual void clear_all()
	{
		for (int i = 0; i < igl_v->data_list.size(); i++)
		{
			igl_v->data_list[i].clear();
		}
	}
	virtual void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int id)
	{
		igl_v->data_list[id].set_mesh(V, F);
	}

	virtual void set_camera_zoom(double zoom)
	{
		igl_v->core().camera_zoom = zoom;
	}

	virtual void set_camera_eye(const Eigen::RowVector3d& p)
	{
		igl_v->core().camera_eye = p.cast<float>();
	}

	virtual void set_camera_center(const Eigen::RowVector3d& c)
	{
		igl_v->core().camera_center = c.cast<float>();
	}
	virtual void set_vertices(const Eigen::MatrixXd& V, int id)
	{
		igl_v->data_list[id].set_vertices(V);
	}
	virtual void compute_normals(int id)
	{
		igl_v->data_list[id].compute_normals();
	}

	virtual void clear(int id)
	{
		igl_v->data_list[id].clear();
	}

	virtual void set_color(const Eigen::RowVector3d & d, int id)
	{
		igl_v->data_list[id].set_colors(d);
	}
	virtual void set_color(const Eigen::MatrixXd& d, int id)
	{
		igl_v->data_list[id].set_colors(d);
	}

	virtual void set_double_sided(bool ds, int id)
	{
		igl_v->data_list[id].double_sided = ds;
	}
	virtual void invert_normals(bool invert_normals, int id)
	{
		igl_v->data_list[id].invert_normals = invert_normals;
	}

	virtual void set_points(const Eigen::MatrixXd& points, int id)
	{
		Eigen::RowVector3d z(0, 0, 0);
		igl_v->data_list[id].set_points(points, z);
	}

	virtual void configure_wireframe_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int id)
	{
		igl_v->data_list[id].clear();
		igl_v->data_list[id].set_mesh(V, F);
		igl_v->data_list[id].show_faces = false;
		igl_v->data_list[id].show_lines = true;
	}

	virtual void set_data(const Eigen::VectorXd &D, int i)
	{
		igl_v->data_list[i].set_data(D);
	}

	virtual void set_data(const Eigen::VectorXd& D,  igl::ColorMapType cmap, int i)
	{
		igl_v->data_list[i].set_data(D, cmap);
	}

	virtual void set_data(const Eigen::VectorXd& D, double min, double max, int num_steps, igl::ColorMapType cmap, int i)
	{
		igl_v->data_list[i].set_data(D, min, max, cmap, num_steps );
	}

	virtual void set_data(const Eigen::VectorXd& D, double min, double max, int i)
	{
		igl_v->data_list[i].set_data(D, min, max);
	}

	virtual void set_points(const Eigen::MatrixXd& points,const  Eigen::RowVector3d& z, int id)
	{
		igl_v->data_list[id].set_points(points, z);
	}

	virtual void set_points(const Eigen::MatrixXd& points,const  Eigen::MatrixXd& color, int id)
	{
		igl_v->data_list[id].set_points(points, color);
	}

	virtual void set_lines(const Eigen::MatrixXd& V,const  Eigen::MatrixXi& E,const  Eigen::MatrixXd& color, int id)
	{
		igl_v->data_list[id].set_edges(V, E, color);
	}

	virtual void set_lines(Eigen::MatrixXd& V,const  Eigen::MatrixXi& E,const  Eigen::RowVector3d& z, int id)
	{
		igl_v->data_list[id].set_edges(V, E, z);
	}

	virtual void set_lines(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, int id)
	{
		Eigen::RowVector3d z(0, 0, 0);
		igl_v->data_list[id].set_edges(V, E, z);
	}
	virtual void set_face_based(bool face_based, int id)
	{
		igl_v->data_list[id].face_based = face_based;
	}
	virtual void set_show_lines(bool show_lines, int id)
	{
		igl_v->data_list[id].show_lines = show_lines;
	}

	virtual bool get_show_lines(int id)
	{
		return igl_v->data_list[id].show_lines;
	}
	virtual void set_show_faces(bool show_faces, int id)
	{
		igl_v->data_list[id].show_faces = show_faces;
	}
	virtual bool get_show_faces(int id)
	{
		return igl_v->data_list[id].show_faces;
	}
	

	virtual void set_data_colormap(const Eigen::VectorXd& d, const Eigen::MatrixXd& cmap, int id)
	{
		igl_v->data_list[id].set_data(d);
		igl_v->data_list[id].set_colormap(cmap);
	}


	virtual void set_texture(const string& tex_png,const  MatrixXd& TC,const  MatrixXi& FTC, int id)
	{
		Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
		bool read = igl::png::readPNG(tex_png, R, G, B, A);
		igl_v->data_list[id].set_colors(Eigen::RowVector3d(1, 1, 1));
		igl_v->data_list[id].set_uv(TC, FTC);
		igl_v->data_list[id].show_texture = true;
		igl_v->data_list[id].set_texture(R, G, B);
		
	}

	//////////////////////// legacy stuff

	virtual void configure_clusters(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters) {};

	virtual void configure_clusters(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters, int id){};

	virtual void configure_deformation_texture(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W){};

	virtual void configure_solid_color_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::RowVector3d& color, int id){};
	virtual void configure_color_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& color, int id){};

	virtual void configure_color_texture(std::string texture_filepath, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
		Eigen::MatrixXd& TC, Eigen::MatrixXi& FTC, int id = 0){};

	virtual void configure_color_texture(std::string texture_filepath, Eigen::MatrixXd& V_coarse, Eigen::MatrixXi& F_coarse, Eigen::MatrixXd& V_fine, Eigen::MatrixXi& F_fine,
		Eigen::MatrixXd& UV_fine, Eigen::MatrixXi& FUV_fine){};

	virtual void configure_matcap(std::string matcap_file, Eigen::MatrixXd& V, Eigen::MatrixXi& F){};

	virtual void set_visible(bool is_visible, int id)
	{
		igl_v->data_list[id].is_visible = is_visible;
	}

	virtual bool get_visible(int id)
	{
		return igl_v->data_list[id].is_visible;
	}
	/*
	Renders mesh  full space mesh V.
	*/
	virtual void render_full(const Eigen::MatrixXd& V, int id) {};


	/*
	Renders mesh with reduced deformation z, and rig parameters sol_p.
	Must do a projection step Bz + Jp to get the final V. Does this step on the CPU
	*/
	virtual void render_reduced_cpu_proj(const Eigen::VectorXd& z,
		const Eigen::VectorXd& p,
		const Eigen::MatrixXd& B,
		const Eigen::SparseMatrix<double>& J, int id) {};

	virtual void render_reduced_gpu_proj(const Eigen::VectorXd& z, 
		const  Eigen::VectorXd& p, int n, int id) {};

	//void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);

	/*
	GUIZMO controls
	*/
	virtual void init_guizmo(bool visible, const Eigen::Matrix4f& A0, std::function<void(const Eigen::Matrix4f& A)>& callback, ImGuizmo::OPERATION op = ImGuizmo::TRANSLATE);

	virtual void set_guizmo_operation(ImGuizmo::OPERATION op)
	{
		guizmo->operation = op;
	};
	ImGuizmo::OPERATION get_guizmo_operation()
	{
		return guizmo->operation;
	}

	virtual void set_animation_max_fps(int max_fps)
	{
		igl_v->core().animation_max_fps = max_fps;
	}
	/*
	configures viewer at index id with viewer parameters sol_p
	*/
	void configure_viewer(fast_cd_viewer_parameters& p, int id)
	{
		invert_normals(p.invert_normals, id);
		set_show_faces(p.show_faces, id);
		set_show_lines(p.show_lines, id);
	}
public:
    igl::opengl::glfw::Viewer v;
	igl::opengl::glfw::Viewer* igl_v;
	igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo; //add guizmo
	igl::opengl::glfw::imgui::ImGuiPlugin* imgui_plugin;

	igl::opengl::glfw::imgui::ImGuiMenu* imgui_menu;
	std::function<bool()> callback;
	int fid, cid; //indices in the viewer.data_list array corresponding to the fine mesh, and the coarse mesh
};
