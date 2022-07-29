#pragma once
#include "SkeletonRig.h"
#include "SkeletonRigFKMouseController.h"
#include "fast_cd_init.h"
#include "ScriptedRigController.h"
/// <summary>
/// Returns the rig controller object corresponding to the rig controller string.
/// </summary>
/// <param name="rig"></param> Rig object. For now, only a HandleRig or a SkeletonRig
/// <param name="rig_controller_str"></param> Can be either "handle", or "skeletonFK" (skeletonIK still needs to be re-exposed)
/// <returns></returns>
RigController* get_rig_controller(std::string rig_controller_str, Rig* rig, FastCDViewer* v);


ScriptedRigController* get_scripted_rig_controller(FastCDInit* init, Rig* rig);