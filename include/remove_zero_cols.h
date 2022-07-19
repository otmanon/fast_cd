#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
/// <summary>
/// Finds 0 columns in A and removes them
/// </summary>
/// <param name="A"></param>
void remove_zero_cols(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, Eigen::VectorXi& I, Eigen::VectorXi& CI);