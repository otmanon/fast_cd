#include "pick_rig_controller.h"
#include "SkeletonRig.h"
#include "HandleRig.h"
#include "SkeletonRigFKMouseController.h"
#include "HandleRigController.h"

RigController * pick_rig_controller(Rig* rig, std::string anim_dir, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
    RigController* controller;
    //picks the default rig controller according to the rig type. 
    if (dynamic_cast<SkeletonRig*>(rig) != nullptr)   //dealing with a skeleton... make a controller
    {
        controller = new SkeletonRigFKMouseController(rig->p0, ((SkeletonRig*)rig)->pI, ((SkeletonRig*)rig)->lengths, viewer, guizmo, anim_dir);
    }
    else
    {
        controller = new HandleRigMouseController(rig->p0, viewer, guizmo, anim_dir);
    }
    return controller;
}