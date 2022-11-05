#include "scale_and_center_geometry.h"

Eigen::MatrixXd scale_and_center_geometry(
   const Eigen::MatrixXd& V, const  double h, const Eigen::RowVector3d c , double& so, Eigen::RowVector3d& to)
{
    Eigen::MatrixXd V2;
    so = 1.0 / (V.col(1).maxCoeff() - V.col(1).minCoeff());
    V2 = V * so;

    to = V2.colwise().mean();
    V2 = V2.rowwise() - to;
    return V2;
}