#include "get_rig.h"
#include "get_rig_type.h"

#include "get_rig_file_format.h"
#include "HandleRig.h"
#include "SkeletonRig.h"
Rig* get_rig(std::string rig_file, Eigen::MatrixXd& V0, Eigen::MatrixXi& T, double radius)
{
	std::string rig_type;
	RIG_TYPE type = get_rig_type(rig_file, rig_type);
    std::string rig_file_format;
    Rig* rig;
    switch (type)
    {
    case RIG_TYPE::LBS_RIG:
        // check if the LBS_rig is in surface format.
        rig_file_format = get_rig_file_format(rig_file);
        if (rig_file_format == "surface")
        {
            rig = new SkeletonRig(rig_file, V0, T, radius);
        }
        else
        {
            rig = new SkeletonRig(rig_file, radius);
        }

        break;
    case RIG_TYPE::AFFINE_RIG:
        rig = new HandleRig(V0, T, radius);
        // as.constraint_controller = new ConstraintController(viewer, guizmo); //null constraints
        break;
    case RIG_TYPE::HANDLE_RIG:
        rig = new HandleRig(rig_file, radius);
        //  as.constraint_controller = new ConstraintController(viewer, guizmo);
         // exit(0);
        break;
    }
}