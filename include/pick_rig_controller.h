#pragma once
#include "rig.h"
#include "RigController.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>

RigController* pick_rig_controller(Rig* rig, std::string anim_dir, igl::opengl::glfw::Viewer * viewer, igl::opengl::glfw::imgui::ImGuizmoWidget * guizmo);