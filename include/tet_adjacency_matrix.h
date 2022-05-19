#pragma once
#include <Eigen/Sparse>

void tet_adjacency_matrix_face_share(const Eigen::MatrixXi& T, Eigen::SparseMatrix<bool, 0, int>& A);


void tet_adjacency_matrix_vert_share(const Eigen::MatrixXi& T, Eigen::SparseMatrix<bool, 0, int>& A);