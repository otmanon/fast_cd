#include "end_effectors_objective_and_gradient.h"
#include "copy_skeleton_at.h"
#include "kinematics_jacobian.h"
#include "transformed_tips.h"

void end_effectors_objective_and_gradient(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    const Eigen::VectorXd& xb0,
    std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f,
    std::function<void(Eigen::VectorXd&)>& proj_z)
{
    f = [&](const Eigen::VectorXd& A)->double
    {
        Skeleton copy = copy_skeleton_at(bone_list, A);
        Eigen::VectorXd xb = transformed_tips(copy, b);
        return 0.5 * (xb - xb0).squaredNorm();
    };
    grad_f = [&](const Eigen::VectorXd& A)->Eigen::VectorXd
    {
        Skeleton copy = copy_skeleton_at(bone_list, A);
        Eigen::VectorXd xb = transformed_tips(copy, b);
        Eigen::VectorXd dEdx = xb - xb0;
        // We want to build dx/dθ where θ are _all_ of the joint angles.
        Eigen::MatrixXd J;
        kinematics_jacobian(copy, b, J);
        return J.transpose() * dEdx;
    };
    proj_z = [&](Eigen::VectorXd& A)
    {
        assert(bone_list.size() * 3 == A.size());
        for (int si = 0; si < bone_list.size(); si++)
        {
            const auto bone = bone_list[si];
            for (int a = 0; a < 3; a++)
            {
                A(si * 3 + a) =
                    std::min(std::max(A(si * 3 + a), bone.xzx_min(a)), bone.xzx_max(a));
            }
        }
    };
}



void end_effectors_kinematic_jacobian(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dxda)
{
    dxda = [&](const Eigen::VectorXd& A)->Eigen::MatrixXd
    {
        Skeleton copy = copy_skeleton_at(bone_list, A);
        // We want to build dx/dθ where θ are _all_ of the joint angles.
        Eigen::MatrixXd J;
        kinematics_jacobian(copy, b, J);
        return J;
    };
}