#pragma once
#include <Eigen/Core>

Eigen::MatrixXd normalize_geometry_unit_height(const Eigen::MatrixXd& V, double& scale)
{
    //ensure mesh is 1 meter tall, centered in the middle

   //V2 = V.rowwise() - V.colwise().mean();
    Eigen::MatrixXd V2;
    scale = 1.0 / (V.col(1).maxCoeff() - V.col(1).minCoeff());
    V2 = V * scale;
    return V2;
}


Eigen::MatrixXd normalize_geometry_unit_height(const Eigen::MatrixXd& V)
{
    //ensure mesh is 1 meter tall, centered in the middle

   //V2 = V.rowwise() - V.colwise().mean();
    Eigen::MatrixXd V2;
    double scale = 1.0 / (V.col(1).maxCoeff() - V.col(1).minCoeff());
    V2 = V * scale;
    return V2;
}

