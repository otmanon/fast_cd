#include "transformed_tips.h"
#include "forward_kinematics.h"
Eigen::VectorXd transformed_tips(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b)
{
    std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> T;
    forward_kinematics(bone_list, T);

    Eigen::Vector3d l;
    Eigen::VectorXd tips = Eigen::VectorXd::Zero(3 * b.size());

    for (int i = 0; i < b.rows(); i++)
    {
        l = Eigen::Vector3d(bone_list[b(i)].length, 0, 0);
        tips.block(3 * i, 0, 3, 1) = T[b(i)] * (bone_list[b(i)].rest_T * l);
    }
    return tips;
}


