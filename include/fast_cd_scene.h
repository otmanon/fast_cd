#pragma once
#include "fast_cd_scene_obj.h"
#include "fast_cd_viewer_custom_shader.h"
#include "fast_cd_viewer_parameters.h";
#include "prolongation.h"
#include <igl/readOBJ.h>
#include <functional>
struct fast_cd_scene
{

	std::vector<fast_cd_scene_obj> scene_objects;
	std::vector<int> object_id;		//object ids organized with scene_objects that finds their data list in the viewer


	fast_cd_viewer_custom_shader viewer;

	std::function<void()> pre_draw_callback;
	fast_cd_scene(std::string vertex_shader, std::string fragment_shader)
	{
		viewer = fast_cd_viewer_custom_shader(vertex_shader, fragment_shader, 16, 16);

		pre_draw_callback = [&]()
		{
			step();
		};
		viewer.set_pre_draw_callback(pre_draw_callback);
	};



	int add_scene_object(fast_cd_scene_obj& obj)
	{
		scene_objects.push_back(obj);
		int id;
		viewer.add_mesh(id);
		object_id.push_back(id);

		viewer.set_mesh(obj.V, obj.F, id);
		viewer.invert_normals(true, 0);
		//where to do we load + attach weights? I think here
		viewer.set_weights(obj.W, obj.Ws, id);//
		return id;
	}


	int add_scene_object(fast_cd_scene_obj& obj, fast_cd_viewer_parameters& p)
	{
		int id;
		viewer.add_mesh(id);
		scene_objects.push_back(obj);
		viewer.configure_viewer(p, id);
		object_id.push_back(id);

		//where to do we load + attach weights? I think here
		if (p.vis_texture)
		{
			MatrixXd V_tex, TC, N;
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
		}
		else
		{
			viewer.set_mesh(obj.V, obj.F, id);
			viewer.set_weights(obj.W, obj.Ws, id);//
		}
		return id;
	}

	// Steps all existing fast_cd_scene_objects once in their animaiton
	void step()
	{
		for (int i = 0; i < scene_objects.size(); i++)
		{
			VectorXd z, p;
			scene_objects[i].step(p, z);
			viewer.set_bone_transforms(p, z, object_id[i]);
			viewer.updateGL(object_id[i]);
		}
	}

	void show(int max_fps)
	{
		viewer.launch(max_fps, true);
	}
};