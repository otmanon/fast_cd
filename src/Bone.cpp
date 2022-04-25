#include "Bone.h"
//Constructor
Bone::Bone(
    int _parent_index,
    int _weight_index,
    Eigen::Affine3d _rest_T,
    double _length, Eigen::VectorXi _wI, Eigen::VectorXd _wV) :
    rest_T(std::move(_rest_T)),
    length(std::move(_length)),
    parent_index(std::move(_parent_index)),
    weight_index(std::move(_weight_index)),
    wI(_wI),
    wV(_wV),
    xzx(0, 0, 0),
    xzx_max(inf, inf, inf),
    xzx_min(-inf, -inf, -inf)
    {
    //need to prerotate T_rest bones by 90 degrees about Z axis 
    //blender makes bones extend in the y axis, not the x axis that we assume here. To get them on the x axis, prerotate by 90 about Z
    //Unsure why these bones then need to be rotated about X again... should look deeper into this if it causes a bigger issue.
    rest_T = Eigen::AngleAxisd(-90 * igl::PI / 180.0, Eigen::Vector3d::UnitX()) * rest_T * Eigen::AngleAxisd(90 * igl::PI / 180.0, Eigen::Vector3d::UnitZ());
    };