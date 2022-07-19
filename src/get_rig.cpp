#include "get_rig.h"
#include "get_rig_type.h"

#include "get_rig_file_format.h"
#include "HandleRig.h"
#include "SkeletonRig.h"
#include "compute_bbw_weights.h"

#include <filesystem>

#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>

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

    return rig;
}



Rig* get_rig(std::string rig_file, Eigen::MatrixXd& V0, Eigen::MatrixXi& T, double radius, std::string weights, std::string weights_type)
{
    namespace fs = std::filesystem;

    std::string rig_type;
    RIG_TYPE type = get_rig_type(rig_file, rig_type);
    std::string rig_file_format;
    Rig* rig;
    
    Eigen::MatrixXd W;
    bool weights_exist = (fs::exists(fs::path(weights)));
    if (weights_exist)
    {
        igl::readDMAT(weights, W);
    }
    if (weights_exist)
    {
        switch (type)
        {
        case RIG_TYPE::LBS_RIG:
            rig = new SkeletonRig(rig_file, W, V0, T, radius);
        
            break;
        case RIG_TYPE::AFFINE_RIG:
            rig = new HandleRig(V0, T, radius);
            break;
        case RIG_TYPE::HANDLE_RIG:
            rig = new HandleRig(rig_file, W, radius);
             // exit(0);
            break;
        }
    }
    else
    {
        Rig* temp_rig;
        switch (type)
        {
        case RIG_TYPE::LBS_RIG:
            // check if the LBS_rig is in surface format.
            rig_file_format = get_rig_file_format(rig_file);
            if (rig_file_format == "surface")
            {
                temp_rig =  new SkeletonRig(rig_file, V0, T, radius);
            }
            else
            {
                temp_rig =  new SkeletonRig(rig_file, radius);
            }
            if (weights_type == "bbw")
            {
                printf("Computing bounded biharmonic weights ...\n");
                W = compute_bbw_weights(temp_rig->p0, V0, T);
                printf("Done!\n");
            }
            else
            {
                W = temp_rig->W;
            }
            igl::writeDMAT(weights, W);
            rig = new SkeletonRig(rig_file, W, V0, T, radius);

            break;
        case RIG_TYPE::AFFINE_RIG:
            rig = new HandleRig(V0, T, radius);
            // as.constraint_controller = new ConstraintController(viewer, guizmo); //null constraints
            break;
        case RIG_TYPE::HANDLE_RIG:
            temp_rig = new HandleRig(rig_file, radius);
            if (rig_type == "bbw")
            {
                printf("Computing bounded biharmonic weights ...\n");
                W = compute_bbw_weights(temp_rig->p0, V0, T);
                printf("Done!\n");
            }
            else
            {
                W = rig->W;
            }
            igl::writeDMAT(weights, W);
            rig = new HandleRig(rig_file, W, radius);
            break;
        }
    }
    return rig;
}