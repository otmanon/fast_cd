#pragma once

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <igl/matlab/matlabinterface.h>

using namespace Eigen;

void matlab_eigs(SparseMatrix<double>& A, SparseMatrix<double>& B, int num_modes, MatrixXd& U, VectorXd& L)
{
    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine, "A", A);
    igl::matlab::mlsetmatrix(&engine, "B", B);
    igl::matlab::mlsetscalar(&engine, "r", num_modes);
    igl::matlab::mleval(&engine, "tic");
    igl::matlab::mleval(&engine, "[EV, S] = eigs(A, B, r,'sm')");
    igl::matlab::mleval(&engine, "EV = real(EV)");
    igl::matlab::mleval(&engine, "S = real(S)");
    igl::matlab::mleval(&engine, "time = toc");
    // igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
    igl::matlab::mlgetmatrix(&engine, "EV", U);
    Eigen::MatrixXd S_mat;
    double t;
    igl::matlab::mlgetmatrix(&engine, "S", S_mat);
    t = igl::matlab::mlgetscalar(&engine, "time");
    printf("Total matlab::eigs decomposition time: %g... \n", t);
    L = S_mat.diagonal(); // Eigen::Map<Eigen::VectorXd>(S_mat.data(), S_mat.rows(), 1);
}