#pragma once

#include "lbs_product_isolate_W.h"
#include "lbs_jacobian.h"
#include "redundant_indices_JV.h"

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <igl/matlab/matlabinterface.h>
#include <igl/repdiag.h>
#include <igl/colon.h>
#include <igl/cat.h>
#include <igl/massmatrix.h>
using namespace Eigen;
using namespace igl;

/*
H -> 3n x 3n full space elastic energy hessian at rest
M -> n x n masssmatrix 
Aeq -> equality constraint matrix c x 3n
num_modes -> num modes to compute
W -> n x w secnodary bone weights
L -> w x 1 secondary bone frequencies

*/
void skinning_modes(const MatrixXd& V, const SparseMatrix<double>& H, const SparseMatrix<double>& M, const SparseMatrix<double>& Aeq, int num_modes, MatrixXd& B_lbs , MatrixXd& W, VectorXd& L)
{
	SparseMatrix<double> S, Sx, Sy, Sz, Z, Q, M2;
	S.resize(H.rows() / 3, H.rows() / 3);
	S.setIdentity();

	Z.resize(S.rows(), S.cols());
	Sx = cat(1, S, cat(1, Z, Z));
	Sy = cat(1, Z, cat(1, S, Z));
	Sz = cat(1, Z, cat(1, Z, S));

	Q = Sx.transpose() * H * Sx + Sy.transpose() * H * Sy + Sz.transpose() * H * Sz;
	

	MatrixXd Weq1;
    lbs_product_isolate_W(V, Aeq.toDense(), Weq1); //confirm that this works somehow

    int dim = V.cols();
    VectorXi I;
    redundant_indices_JV(Aeq.rows() / (dim * (dim + 1)), V.cols(), I);
    MatrixXi Im = MatrixXi(I.rows(), 1);
    Im = I;
    SparseMatrix<double> Weq = Weq1.sparseView();

    Engine* engine;
    matlab::mlinit(&engine);
    matlab::mlsetmatrix(&engine, "A", Q);
    matlab::mlsetmatrix(&engine, "B", M);
    matlab::mlsetmatrix(&engine, "Aeq", Weq);
    matlab::mlsetmatrix(&engine, "redI", Im);
    //remove redundant rows fiiirst
 //   matlab::mleval(&engine, "[~, p] = rref(Aeq')");
    matlab::mleval(&engine, "Aeq = Aeq(redI,  :)");


    matlab::mlsetscalar(&engine, "r", num_modes);
    matlab::mleval(&engine, "tic");
    matlab::mleval(&engine, "A = [A Aeq'; Aeq zeros(size(Aeq, 1))];");
    matlab::mleval(&engine, "M = blkdiag(B, zeros(size(Aeq, 1)))");
    matlab::mleval(&engine, "[EV, S] = eigs(A + speye(size(A))*1e-8, M, r,'sm')");
    matlab::mleval(&engine, "EV = real(EV)");
    matlab::mleval(&engine, "S = real(S)");
    matlab::mleval(&engine, "time = toc");
    // igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
    MatrixXd W2;
    igl::matlab::mlgetmatrix(&engine, "EV", W2);
    Eigen::MatrixXd S_mat;
    double t;
    igl::matlab::mlgetmatrix(&engine, "S", S_mat);
    t = igl::matlab::mlgetscalar(&engine, "time");
    printf("Total matlab::eigs decomposition time: %g... \n", t);
    L = S_mat.diagonal(); // Eigen::Map<Eigen::VectorXd>(S_mat.data(), S_mat.rows(), 1);

    W = W2.topRows(V.rows());
    SparseMatrix<double> B_test;
    lbs_jacobian(V, W, B_test);
    B_lbs = B_test.toDense();
    //Eigen::MatrixXd err = Aeq * B_test;
    //std::cout << "error " << err <<  std::endl;
   // assert((Aeq * B_test).norm() < 1e-6);
}