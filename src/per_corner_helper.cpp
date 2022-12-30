#include "per_corner_helper.h"

void per_corner_helper(const Eigen::MatrixXd& X, const Eigen::MatrixXi& F,
    igl::opengl::MeshGL::RowMatrixXf& X_vbo)
{
    X_vbo.resize(F.rows() * 3, X.cols());
    for (unsigned i = 0; i < F.rows(); ++i)
        for (unsigned j = 0; j < 3; ++j)
            X_vbo.row(i * 3 + j) = X.row(F(i, j)).cast<float>();

}