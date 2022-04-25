#include "euler_angles_to_transform.h"
#include <igl/PI.h>
Eigen::Affine3d euler_angles_to_transform(
    const Eigen::Vector3d& xzx)
{
     /////////////////////////////////////////////////////////////////////////////
  // This is faster than doing it in the condense way
     
     Eigen::Affine3d R;
     {
         R.matrix() <<
             1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
         Eigen::Vector3d rad_a = EIGEN_PI / 180 * xzx;
         Eigen::Matrix3d rx1, rz, rx2;
         rx1 << 1, 0, 0,
             0, cos(rad_a[0]), -sin(rad_a[0]),
             0, sin(rad_a[0]), cos(rad_a[0]);

         rx2 << 1, 0, 0,
             0, cos(rad_a[2]), -sin(rad_a[2]),
             0, sin(rad_a[2]), cos(rad_a[2]);

         rz << cos(rad_a[1]), -sin(rad_a[1]), 0,
             sin(rad_a[1]), cos(rad_a[1]), 0,
             0, 0, 1;

         R.matrix().block(0, 0, 3, 3) = rx1 * rz * rx2;
     }

     return R;

    /////////////////////////////////////////////////////////////////////////////
}
