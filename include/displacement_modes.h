#pragma once
#include <augment_with_linear_constraints.h>

#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <igl/get_seconds.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/cat.h>
#include <igl/blkdiag.h>
using namespace Eigen;

void displacement_modes(const SparseMatrix<double>& H, const SparseMatrix<double>& M, const int num_modes,
	MatrixXd& B, VectorXd& V)
{
    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine, "A", H);
    igl::matlab::mlsetmatrix(&engine, "B", M);
    igl::matlab::mlsetscalar(&engine, "r", num_modes);
    igl::matlab::mleval(&engine, "tic");
    igl::matlab::mleval(&engine, "[EV, S] = eigs(A, B, r,'sm')");
    igl::matlab::mleval(&engine, "EV = real(EV)");
    igl::matlab::mleval(&engine, "S = real(S)");
    igl::matlab::mleval(&engine, "time = toc");
    // igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
    igl::matlab::mlgetmatrix(&engine, "EV", B);
    Eigen::MatrixXd S_mat;
    double t;
    igl::matlab::mlgetmatrix(&engine, "S", S_mat);
    t = igl::matlab::mlgetscalar(&engine, "time");
    printf("Total matlab::eigs decomposition time: %g... \n", t);
    V = S_mat.diagonal(); // Eigen::Map<Eigen::VectorXd>(S_mat.data(), S_mat.rows(), 1);
}

/// <summary>
/// Computes normal linear displacement modes for an elastic energy with hessian H, and constraints described by J
/// </summary>
/// <param name="H">d*n x d*n elastic energy hessian </param>
/// <param name="M"> d*n x d*n massmatrix </param>
/// <param name="J"> #c x d*n constraint matrix</param>
/// <param name="num_modes"> number of modes </param>
/// <param name=""> d*n x num_modes subspace matrix (eigenvectors) </param>
/// <param name=""> num_modes x 1 subspace modal frequencies (eigenvalues) </param>
void displacement_modes(const SparseMatrix<double>& H, const SparseMatrix<double>& M, const SparseMatrix<double>& J, const int num_modes, 
	MatrixXd& B, VectorXd& V)
{
    Eigen::SparseMatrix<double> Q, MZ;
    augment_with_linear_constraints(H, J, Q);

    Eigen::SparseMatrix<double> Z = Eigen::SparseMatrix<double>(J.rows(), J.rows());
    igl::blkdiag({ M, Z }, MZ);

    Eigen::MatrixXd B2;
    displacement_modes(Q, MZ, num_modes, B2, V);

    B = B2.block(0, 0, H.rows(), B2.cols());
}