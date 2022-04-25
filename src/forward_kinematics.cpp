#include "forward_kinematics.h"
#include "euler_angles_to_transform.h"
/*
Applies affine transformation to each bone in order to get the world space positions
*/
void forward_kinematics(const Skeleton& bone_list, std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T)
{

    T.resize(bone_list.size(), Eigen::Affine3d::Identity());
    Eigen::Affine3d local_rot, T_hat_i, T_hat_inv_i, T_i, T_p;

    std::vector<bool> computed;
    computed.resize(T.size(), false);

    std::function<void(int) > get_transform = [&](int b)
    {

        if (!computed[b])
        {
            if (bone_list[b].parent_index < 0 && !computed[b])
            {
                //sets parent bone to no rotation, no translation
                T[b].setIdentity();
            }
            else
            {
                get_transform(bone_list[b].parent_index);
                local_rot = euler_angles_to_transform(bone_list[b].xzx);
                T_hat_i = bone_list[b].rest_T;
                T_hat_inv_i = bone_list[b].rest_T.inverse();

                T[b] = T[bone_list[b].parent_index] * (T_hat_i * (local_rot * (T_hat_inv_i)));
            }
            computed[b] = true;
        }
    };

    for (int i = 0; i < bone_list.size(); i++)
    {
        get_transform(i);
    }

}

