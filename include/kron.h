#pragma once
#include <Eigen/Core>

/// <summary>
/// Performs C = kron(A, B), where C_ij = [A B_ij]
/// </summary>
/// <typeparam name="DerivedA"></typeparam>
/// <typeparam name="DerivedB"></typeparam>
/// <typeparam name="DerivedC"></typeparam>
/// <param name="A"></param>
/// <param name="B"></param>
/// <param name="C"></param>

using namespace Eigen;
void kron(const MatrixXd& A, const MatrixXd& B, MatrixXd& C);
