#pragma once
#include "fast_cd_scene_obj.h"
#include "fast_cd_viewer_custom_shader.h"
#include "fast_cd_viewer_parameters.h";
#include "prolongation.h"
#include "ScalarRecorder.h"

#include "VectorRecorder.h"


#include <igl/writeOBJ.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/readOBJ.h>
#include <igl/ray_mesh_intersect.h>
#include <filesystem>
#include <functional>
struct fast_cd_scene
{

	std::vector<fast_cd_scene_obj> scene_objects;

	//std::vector<MeshRecorder> mesh_recorders;
	std::vector<VectorRecorder> z_recorders, p_recorders;



	bool record_meshes;
	std::string record_dir;

	std::vector<int> object_id;		//object ids organized with scene_objects that finds their data list in the viewer

	fast_cd_viewer_custom_shader viewer;

	std::function<void()> pre_draw_callback;
	std::function<void()> menu_callback;
	std::function<void(int button, int modifie)> mouse_down_callback;
	std::function<bool(unsigned int key, int modifier)> key_pressed_callback;
	std::function<void(const Matrix4f& )> guizmo_callback;
	//std::function<void()> mouse_up_callback;
	//std::function<void()> mouse_pressed_callback;

	ScalarRecorder totalt;
	int num_v, num_t, num_obj;

	double timestep;
	bool do_cd; //if this is on then overrides all objects do_cd

	int max_anim_length;


	int controllable_object_id;
	VectorXd p_frozen;
	Matrix4f A0;
	fast_cd_scene(std::string vertex_shader, std::string fragment_shader, int num_primary_bones = 16, int num_secondary_bones = 16)
	{
		viewer = fast_cd_viewer_custom_shader(vertex_shader, fragment_shader, num_primary_bones, num_secondary_bones);

		pre_draw_callback = [&]()
		{
			step();
		};

		key_pressed_callback = [&](unsigned int key, int modifier)->bool
		{
			if (key == 'g')
			{
				if (ImGuizmo::TRANSLATE == viewer.get_guizmo_operation())
				{
					viewer.set_guizmo_operation(ImGuizmo::ROTATE);
				}
				else if (ImGuizmo::ROTATE == viewer.get_guizmo_operation())
				{
					viewer.set_guizmo_operation(ImGuizmo::TRANSLATE);
				}
			}

			return false;
		};

		viewer.set_key_pressed_callback(key_pressed_callback);
		A0.setIdentity();
		guizmo_callback = [&](const Matrix4f& T)
		{
			MatrixXd A = (A0.inverse() * T).topRows(3).cast<double>();

			VectorXd p_tmp = p_frozen;
			transform_rig_parameters(p_tmp, A);
			scene_objects[controllable_object_id].p_controller = p_tmp;
		};
		A0 = Matrix4f::Identity();
		viewer.init_guizmo(false, A0, guizmo_callback, ImGuizmo::TRANSLATE);

		mouse_down_callback = [&](int button, int modifier)
		{
			Eigen::RowVector3f last_mouse = Eigen::RowVector3f(
				viewer.igl_v->current_mouse_x, viewer.igl_v->core().viewport(3) - viewer.igl_v->current_mouse_y, 0);
			if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Left)
				printf("Left Click\n");

			else if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Right)
				printf("Right Click \n");//

			printf("button %d \n", button);
			if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Right)  //if right click, we are placing handles
			{
				
				int  old_controllable_object_id = controllable_object_id;
				controllable_object_id = pick_scene_object(last_mouse);
				if (controllable_object_id >= 0)
				{

					if (controllable_object_id != old_controllable_object_id && old_controllable_object_id >= 0)
					{
						printf("old controllable object id %d \n", old_controllable_object_id);
						scene_objects[old_controllable_object_id].controlled = false;
						MatrixXd A = (A0.inverse() * viewer.guizmo->T).topRows(3).cast<double>();
						scene_objects[old_controllable_object_id].transform_animation(A);
					}

					scene_objects[controllable_object_id].controlled = true;
					//initialize  A0
					A0 = Matrix4f::Identity();

					//calculate the  center of mass. 
					VectorXd u =
						scene_objects[controllable_object_id].sim.params->B * scene_objects[controllable_object_id].st.z_curr +
						scene_objects[controllable_object_id].sim.params->J * scene_objects[controllable_object_id].st.p_curr;

					MatrixXd U = Map<MatrixXd>(u.data(), u.rows() / 3, 3);
					if (scene_objects[controllable_object_id].do_texture)
						U = scene_objects[controllable_object_id].P* U;
					RowVector3d center = U.colwise().mean();
					A0.block(0, 3, 3, 1) = center.transpose().cast<float>();

					p_frozen = scene_objects[controllable_object_id].st.p_curr;
					//TODO what is p_frozen;
					scene_objects[controllable_object_id].p_controller = p_frozen;;
					set_guizmo_visible(true);
					viewer.guizmo->T = A0;
					printf("id %i, \n", controllable_object_id);
				}
				else if (old_controllable_object_id >= 0)
					controllable_object_id = old_controllable_object_id; 
			}
		};

		viewer.set_pre_draw_callback(pre_draw_callback);

		viewer.set_mouse_down_callback(mouse_down_callback);


		double full_t = 0, sim_t = 0;

		num_obj = 0;
		num_v = 0;
		num_t = 0;

		do_cd = true;
		record_meshes = false;
		timestep = 0;
		max_anim_length = 0;
		menu_callback = [&]()
		{
			ImGui::Text(" #Vertices : %i \n #Tets : %i \n  #Obj : %i \n FPS : %g \n", num_v,
				num_t, num_obj, 1.0/totalt.mean(24));
		};
		viewer.set_menu_callback(menu_callback);

		controllable_object_id = -1;
		p_frozen = VectorXd(0);
		A0.setZero();
	};
	


	int pick_scene_object(Eigen::RowVector3f& last_mouse)
	{
		// Find closest point on mesh to mouse position
		int hit_id = -1;
		for (int i = 0; i < scene_objects.size(); i++)
		{
			int fid;
			Eigen::Vector3f bary;
			VectorXd u =
				scene_objects[i].sim.params->B * scene_objects[i].st.z_curr +
				scene_objects[i].sim.params->J * scene_objects[i].st.p_curr;

			MatrixXd U =  Map<MatrixXd>(u.data(), u.rows() / 3, 3);
			
			if (scene_objects[i].do_texture)
				U = scene_objects[i].P* U;
			MatrixXi F = viewer.igl_v->data_list[object_id[i]].F;
			if (igl::unproject_onto_mesh(
				last_mouse.head(2),
				viewer.igl_v->core().view,
				viewer.igl_v->core().proj,
				viewer.igl_v->core().viewport,
				U, F,
				fid, bary))
			{
				hit_id = i;
				break;
			}
		}

		return hit_id;
	}

	void set_guizmo_visible(bool visibility)
	{
		viewer.guizmo->visible = visibility;
	}

	void set_do_cd(bool do_cd)
	{
		this->do_cd = do_cd;
	}

	void set_background_color(const Vector3d& color)
	{
		viewer.igl_v->core().background_color.topRows(3) = color.cast<float>();
	}


	int add_scene_object(fast_cd_scene_obj& obj)
	{
		scene_objects.push_back(obj);
		int id;
		viewer.add_mesh(id);
		object_id.push_back(id);

		viewer.set_mesh(obj.V, obj.F, id);
		viewer.invert_normals(true, id);
		//where to do we load + attach weights? I think here
		viewer.set_weights(obj.W, obj.Ws, id);//

		num_obj += 1;
		num_v += obj.V.rows();
		num_t += obj.T.rows();
		
		VectorRecorder z_rec(obj.sim.params->B.cols()), p_rec(obj.sim.params->J.cols());
		z_recorders.push_back(z_rec);
		p_recorders.push_back(p_rec);
		return id;
	}


	int add_scene_object(fast_cd_scene_obj& obj, fast_cd_viewer_parameters& p)
	{
		int id;
		viewer.add_mesh(id);
		viewer.configure_viewer(p, id);
		//where to do we load + attach weights? I think here
		if (p.vis_texture)
		{
			MatrixXd TC, N, V_tex;
			MatrixXi F_tex, FTC, FN;
			readOBJ(p.texture_obj, V_tex, TC, N, F_tex, FTC, FN);
			V_tex = (V_tex * p.so).rowwise() - p.to;
			viewer.set_mesh(V_tex, F_tex, id);
			viewer.set_face_based(false, id);
			viewer.set_texture(p.texture_png, TC, FTC, id);
			
			SparseMatrix<double> Pr;
			prolongation(V_tex, obj.V,  obj.T, Pr);
			MatrixXd W_tex = Pr * obj.W;
			MatrixXd Ws_tex = Pr * obj.Ws;
			viewer.set_weights(W_tex, Ws_tex, id);
			
			obj.do_texture = true;
			obj.P = Pr;
			obj.V_tex = V_tex;
			obj.N = N;
			obj.TC = TC;
			obj.F_tex = F_tex;
			obj.FN = FN;
			obj.FTC = FTC;

		}
		else
		{
			viewer.set_mesh(obj.V, obj.F, id);
			viewer.set_weights(obj.W, obj.Ws, id);//
		}

		scene_objects.push_back(obj);
		object_id.push_back(id);

		num_obj += 1;
		num_v += obj.V.rows();
		num_t += obj.T.rows();
		VectorRecorder z_rec(obj.sim.params->B.cols()), p_rec(obj.sim.params->J.cols());
		z_recorders.push_back(z_rec);
		p_recorders.push_back(p_rec);
		max_anim_length = max(max_anim_length, obj.anim_length);

		std::cout << obj.do_texture << ": do_texture" << std::endl;
		return id;
	}

	// Steps all existing fast_cd_scene_objects once in their animaiton
	void step()
	{
		double tstart_sim = igl::get_seconds();
		for (int i = 0; i < scene_objects.size(); i++)
		{
			VectorXd z, p;
			scene_objects[i].step(p, z);
			if (!do_cd)
				z.setZero();
			viewer.set_bone_transforms(p, z, object_id[i]);
			viewer.updateGL(object_id[i]);
			
			if (record_meshes)
			{
				z_recorders[i].record_frame(z);
				p_recorders[i].record_frame(p);
			}
		}
		double tend_sim = igl::get_seconds();
		totalt.record_frame(tend_sim - tstart_sim);
		if (timestep == max_anim_length && record_meshes)
		{
			save_mesh_screenshots();
		}
		timestep += 1;
	}

	void save_mesh_screenshots()
	{
		if (record_meshes)
		{
		for (int i = 0; i < scene_objects.size(); i++)
		{
				VectorRecorder z_rec, p_rec;
				z_rec = z_recorders[i];
				p_rec = p_recorders[i];
				MatrixXi T = scene_objects[i].T;
				MatrixXi F = scene_objects[i].F;
				for (int j = 0; j < z_rec.P.cols(); j++)
				{
					VectorXd z, p;
					z = z_rec.P.col(j);
					p = p_rec.P.col(j);
					
					fast_cd_scene_obj& obj = scene_objects[i];
					
					VectorXd v = (obj.sim.params->B * z +
						obj.sim.params->J * p);
					MatrixXd V = Map<MatrixXd>(v.data(), v.rows() / 3, 3);
			
					std::string dir = record_dir + "/mesh_recordings/object_"
						+ obj.name + "/";
					
					if (obj.do_texture)
					{
						printf("Saving with texture !\n");
						MatrixXd V_tex =  obj.P * V;
						igl::writeOBJ(dir + "/" + std::to_string(j) + ".obj",
						V_tex, obj.F_tex , obj.N, obj.FN, obj.TC, obj.FTC);

					}
					else
					{
						printf("Saving without texture !\n");
						igl::writeOBJ(dir + "/" + std::to_string(j) + ".obj",
							V, F);
					}
					
					namespace fs = std::filesystem;
					if (!fs::exists(fs::path(dir)))
						fs::create_directories(fs::path(dir));

				}
			}
		}
	}

	void show(int max_fps)
	{
		viewer.launch(max_fps, true);
	}

	void set_record(bool record_meshes, std::string dir)
	{
		this->record_meshes = true;
		this->record_dir = dir;
	}
};