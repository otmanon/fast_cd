#include "fast_cd_viewer.h"
#include "launch_viewer_custom_shader.h"
#include "rainbow_cmap.h"

#include <igl/LinSpaced.h>
#include <igl/opengl/gl.h>
#include <igl/png/readPNG.h>
fast_cd_viewer::fast_cd_viewer()
{

	igl_v = & this->v;

	//while (igl_v->data_list.size() < 3) igl_v->append_mesh();

    igl_v->core().background_color.setOnes();
	igl_v->core().animation_max_fps = 60;
	igl_v->core().is_animating = true;

  
    igl_v->callback_key_pressed = [&](igl::opengl::glfw::Viewer& viewer, unsigned int unicode_key, int modifiers)->bool
    {
        return this->default_key_pressed_callback(viewer, unicode_key, modifiers, 0);
    };

    imgui_plugin = new igl::opengl::glfw::imgui::ImGuiPlugin();
    guizmo = new igl::opengl::glfw::imgui::ImGuizmoWidget();
 //   imgui_menu = new igl::opengl::glfw::imgui::ImGuiMenu();


    igl_v->plugins.push_back(imgui_plugin);
    // push back menu here

    imgui_plugin->widgets.push_back(guizmo);
  //  imgui_plugin->widgets.push_back(imgui_menu);
    guizmo->visible = false;

    imgui_menu = new igl::opengl::glfw::imgui::ImGuiMenu();
    imgui_plugin->widgets.push_back(imgui_menu);
    imgui_menu->callback_draw_viewer_menu = [&]()
    {    };
  //  imgui_menu->draw_custom_window = ()[]
}




bool fast_cd_viewer::default_key_pressed_callback(igl::opengl::glfw::Viewer& viewer, unsigned int unicode_key, int modifiers, int id)
{
    
    switch(unicode_key)
    {
      case 'A':
      case 'a':
      {
        igl_v->core().is_animating = !igl_v->core().is_animating;
        return true;
      }
      case 'D':
      case 'd':
      {
        igl_v->data_list[id].double_sided = !igl_v->data_list[id].double_sided;
        return true;
      }
      case 'F':
      case 'f':
      {
          igl_v->data_list[id].set_face_based(!igl_v->data_list[id].face_based);
          return true;
      }
      case 'I':
      case 'i':
      {
          igl_v->data_list[id].dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
          igl_v->data_list[id].invert_normals = !igl_v->data_list[id].invert_normals;
          return true;
      }
      case 'L':
      case 'l':
      {
          igl_v->core().toggle(igl_v->data_list[id].show_lines);
          return true;
      }
      case 'O':
      case 'o':
      {
          igl_v->core().orthographic = !igl_v->core().orthographic;
          return true;
      }
      case 'T':
      case 't':
      {
          igl_v->core().toggle(igl_v->data_list[id].show_faces);
        return true;
      }
    
    }
    return false;
  
}

fast_cd_viewer::fast_cd_viewer(igl::opengl::glfw::Viewer* igl_v, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
    this->igl_v = igl_v;
    while (igl_v->data_list.size() < 3) igl_v->append_mesh();
    igl_v->core().animation_max_fps = 240;
    igl_v->core().is_animating = true;
    fid = 1;
    cid = 0;
    this->guizmo = guizmo;
}


void fast_cd_viewer::init_guizmo(bool visible, const Eigen::Matrix4f& A0,  std::function<void(const Eigen::Matrix4f& A)>& callback, ImGuizmo::OPERATION op)
{
    guizmo->operation = op;
    guizmo->visible = visible;
    guizmo->T = A0;
    guizmo->callback = callback;

}



fast_cd_viewer::fast_cd_viewer(igl::opengl::glfw::Viewer* igl_v)
{
    this->igl_v = igl_v;
	while (igl_v->data_list.size() < 3) igl_v->append_mesh();
    igl_v->core().animation_max_fps = 240;
    igl_v->core().is_animating = true;
    fid = 1;
    cid = 0;

    guizmo = new igl::opengl::glfw::imgui::ImGuizmoWidget();
    imgui_plugin = new igl::opengl::glfw::imgui::ImGuiPlugin();
    igl_v->plugins.push_back(imgui_plugin);
    // push back menu here
    imgui_plugin->widgets.push_back(guizmo);

}

void fast_cd_viewer::launch()
{
    // std::cout << "reach here" << std::endl;
	// launch_viewer_custom_shader(*igl_v, true, false, "fast CD app", 1920, 1080);
    
    igl_v->launch(true, false, "fast cd app", 1920, 1080);
}



void fast_cd_viewer::set_pre_draw_callback(std::function<void()>& callback)
{
	igl_v->callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
	{
		callback();
		return false;
	};
}

