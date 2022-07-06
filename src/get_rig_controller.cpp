#include "get_rig_controller.h"

RigController * get_rig_controller(std::string rig_controller_str, Rig* rig, FastCDViewer* v)
{
	if (rig_controller_str == "handle")
	{
		return new HandleRigMouseController(rig->p0, v);
	}
	else if (rig_controller_str == "skeletonFK" && dynamic_cast<SkeletonRig*>(rig) != nullptr)
	{
		return new SkeletonRigFKMouseController(((SkeletonRig*)rig)->p0, ((SkeletonRig*)rig)->pI, ((SkeletonRig*)rig)->lengths, v);
	}
	else
	{
		return new HandleRigMouseController(rig->p0, v);
	}
}