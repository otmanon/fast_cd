#pragma once
#include "fast_cd_scene_obj.h"
#include "fast_cd_viewer_custom_shader.h"
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

	void add_scene_object(fast_cd_scene_obj& obj)
	{
		scene_objects.push_back(obj);
		int id;
		viewer.add_mesh(id);
		object_id.push_back(id);

		viewer.set_mesh(obj.V, obj.F, id);
		viewer.invert_normals(true, 0);
		//where to do we load + attach weights? I think here
		viewer.set_weights(obj.W, obj.Ws, id);//
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