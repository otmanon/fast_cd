#pragma once
#include <eigen/core>
#include <eigen/sparse>

using namespace Eigen;
/*
Given a vector of values, build a sparse diagonal matrix S from that vector such that S[i, i] = diag[i]
This is mainly implemented because Eigen weirdly cant convert from a diagonal matrix to a sparse matrix efficiently.
*/
SparseMatrix<double> sparse_diag(VectorXd& diag);