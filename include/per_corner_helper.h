#pragma once
#include <igl/opengl/MeshGL.h>

//if face based
  // Input:
//   X  #V by dim quantity
// Output:
//   X_vbo  #F*3 by dim scattering per corner
void per_corner_helper(const Eigen::MatrixXd& X, const Eigen::MatrixXi& F,
    igl::opengl::MeshGL::RowMatrixXf& X_vbo);