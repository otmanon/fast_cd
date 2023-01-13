#pragma once
#include "fast_cd_scene_obj.h"
#include "fast_cd_viewer_custom_shader.h"
#include "fast_cd_viewer_parameters.h";
#include "prolongation.h"
#include "ScalarRecorder.h"

#include "VectorRecorder.h"


#include <igl/writeOBJ.h>

#include <igl/readOBJ.h>
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
	ScalarRecorder totalt;
	int num_v, num_t, num_obj;

	double timestep;
	bool do_cd; //if this is on then overrides all objects do_cd

	int max_anim_length;
	fast_cd_scene(std::string vertex_shader, std::string fragment_shader, int num_primary_bones = 16, int num_secondary_bones = 16)
	{
		viewer = fast_cd_viewer_custom_shader(vertex_shader, fragment_shader, num_primary_bones, num_secondary_bones);

		pre_draw_callback = [&]()
		{
			step();
		};

		viewer.set_pre_draw_callback(pre_draw_callback);
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

	};
	
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