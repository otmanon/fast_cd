#include "compute_modes_matlab.h"
#include <iostream>
#include <cassert>
#include <igl/sort.h>
#include <igl/get_seconds.h>
#include <igl/matlab/matlabinterface.h>


void compute_modes_matlab(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, int r, Eigen::MatrixXd& U, Eigen::VectorXd& S)
{
    // Matlab instance
    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine, "A", A);
    igl::matlab::mlsetmatrix(&engine, "B", B);
    igl::matlab::mlsetscalar(&engine, "r", r);
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
    S = Eigen::Map<Eigen::VectorXd>(S_mat.data(), S_mat.rows(), 1);
  
}
