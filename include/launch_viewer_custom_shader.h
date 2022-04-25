#pragma once
#include <igl/opengl/glfw/Viewer.h>
void launch_viewer_custom_shader(igl::opengl::glfw::Viewer& v,
	bool resizeable=true, bool fullscreen=false, std::string name="app", int width = 1920, int height = 1080);