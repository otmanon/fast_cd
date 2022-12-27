#include "Eigen/Core"
#pragma once
using namespace std;
using namespace Eigen;
struct fast_cd_viewer_parameters
{
	
	bool invert_normals = true;
	bool show_faces = true;
	bool show_lines = false;
	
	bool vis_texture = false;
	string texture_obj; //fill this up if we set vis_texture=true
	string texture_png;
	double so;
	RowVector3d to;

	void set_texture(string texture_obj, string texture_png)
	{
		vis_texture = true;
		this->texture_obj = texture_obj;
		this->texture_png = texture_png;
		this->so = 1.0;
		this->to = RowVector3d(0, 0, 0);
	}

	void set_texture(string texture_obj, string texture_png, double so, RowVector3d to)
	{
		vis_texture = true;
		invert_normals = false;
		this->texture_obj = texture_obj;
		this->texture_png = texture_png;
		this->so = so;
		this->to = to;
	}
};