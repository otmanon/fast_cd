#include "copy_skeleton_at.h"
Skeleton copy_skeleton_at(
    const Skeleton& bone_list,
    const Eigen::VectorXd& A)
{

    Skeleton copy = bone_list;
    for (int b = 0; b < bone_list.size(); b++)
    {
        copy[b].xzx = A.block(3 * b, 0, 3, 1);
    }
    return copy;
}