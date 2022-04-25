#include "kinematics_jacobian.h"
#include "transformed_tips.h"
void kinematics_jacobian(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    Eigen::MatrixXd& J)
{
    // Copy current configuration
    //shouldnt need to copy this for every change you make... just perturb it and then delete it again
    Skeleton copy = bone_list;

    J = Eigen::MatrixXd::Zero(b.size() * 3, bone_list.size() * 3);
    // loop over each bone
    for (int si = 0; si < bone_list.size(); si++)
    {
        // loop over each angle
        for (int a = 0; a < 3; a++)
        {
            // central differences requires one forward, one backward
            for (int dir = -1; dir <= 1; dir += 2)
            {

                // small change
                const double epsilon = 1e-7;
                copy[si].xzx(a) += double(dir) * epsilon;
                Eigen::VectorXd xb = transformed_tips(copy, b);
                for (int bi = 0; bi < b.size(); bi++)
                {
                    J.block(bi * 3, si * 3 + a, 3, 1) +=
                        double(dir) * xb.block(bi * 3, 0, 3, 1) / (2.0 * epsilon);
                }
                //now that you're done with it
                //just undo the tiny change to current configuration... 
                //instead of copying allll the bones over again
                copy[si].xzx(a) -= double(dir) * epsilon;
            }
        }
    }
}
