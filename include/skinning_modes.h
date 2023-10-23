#pragma once

#include "lbs_product_isolate_W.h"
#include "lbs_jacobian.h"
#include "redundant_indices_JV.h"
#include "write_sparse_ijv_DMAT.h"

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <igl/matlab/matlabinterface.h>
#include <igl/repdiag.h>
#include <igl/colon.h>
#include <igl/cat.h>
#include <igl/massmatrix.h>
#include <igl/slice.h>

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
void skinning_modes(const MatrixXd& V, const SparseMatrix<double>& H, const SparseMatrix<double>& Mw, const SparseMatrix<double>& J, int num_modes, MatrixXd& B_lbs, MatrixXd& W, VectorXd& L, bool debug = false, string output_dir = "")
{

    //TODO:: should make a weight space hessian method. The input to this function should always be Hw
    
    SparseMatrix<double> K; // proxy hessian
    if (J.rows() == 0)
    {
        SparseMatrix<double> Id(H.rows(), H.cols());
        Id.setIdentity();
        K = H + 1e-8 * Id; // a little regularization if we have no null space!
    }
    else
        K = H;
    //weight space hessian
    //TODO: make weight_space_hessian(H, Hw);
    SparseMatrix<double> S, Sx, Sy, Sz, Z, Hw;
    S.resize(H.rows() / 3, H.rows() / 3);
    S.setIdentity();

    Z.resize(S.rows(), S.cols());
    Sx = cat(1, S, cat(1, Z, Z));
    Sy = cat(1, Z, cat(1, S, Z));
    Sz = cat(1, Z, cat(1, Z, S));
    Hw = Sx.transpose() * K * Sx + Sy.transpose() * K * Sy + Sz.transpose() * K * Sz;

    // form weight space constraint.
    MatrixXd Jw_redundant, Jwd;
    MatrixXd Jd = J.toDense();
    lbs_product_isolate_W(V, Jd, Jw_redundant); 

    VectorXi I; //redundant indicies
    int dim = V.cols();
    int num_b = J.rows() / (dim * (dim + 1));
  //  redundant_indices_JV(num_b, dim, I);
  //  igl::slice(Jw_redundant, I, 1, Jwd);
    MatrixXi Im = MatrixXi(I.rows(), 1); //need to convert to matrix to pass to matlab... annoying
    Im = I;

  
    //igl::slice(Weq, Im, 2, Jw);
    SparseMatrix<double> Jw= Jw_redundant.sparseView();

    Engine* engine;
    matlab::mlinit(&engine);
    matlab::mlsetmatrix(&engine, "Hw", Hw);
    matlab::mlsetmatrix(&engine, "Mw", Mw);
    matlab::mlsetmatrix(&engine, "Jw", Jw);
   // matlab::mlsetmatrix(&engine, "redI", Im);
    //remove redundant rows fiiirst
    matlab::mleval(&engine, "rows_with_nonzeros = find(~all(Jw == 0, 2));");
    matlab::mleval(&engine, "Jw = Jw(rows_with_nonzeros, :);");

    matlab::mleval(&engine, "[~, p] = rref(full(Jw'))");
    matlab::mleval(&engine, "Jw = Jw(p,  :)");
    matlab::mlsetscalar(&engine, "r", num_modes);
    matlab::mleval(&engine, "tic");
    matlab::mleval(&engine, "H = [Hw Jw'; Jw zeros(size(Jw, 1))];");
    matlab::mleval(&engine, "M = blkdiag(Mw, zeros(size(Jw, 1)))");
    matlab::mleval(&engine, "[EV, S] = eigs(H, M, r,'sm')");
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



    if (debug)
    {
        if (output_dir == "")
        {
            printf("output_dir is empty, have no place to store debug info to! \n");
        }
        else
        {
            write_sparse_ijv_DMAT(output_dir + "/Mw.DMAT", Mw);
            write_sparse_ijv_DMAT(output_dir + "/H.DMAT", H);
            write_sparse_ijv_DMAT(output_dir + "/Aeq.DMAT", J);
            write_sparse_ijv_DMAT(output_dir + "/Hw.DMAT", Hw);
            write_sparse_ijv_DMAT(output_dir + "/Jw.DMAT", Jw);
            write_sparse_ijv_DMAT(output_dir + "/J.DMAT", J);
            igl::writeDMAT(output_dir + "/rI.DMAT", I);
            igl::writeDMAT(output_dir + "/W_cpp.DMAT", W);
            igl::writeDMAT(output_dir + "/L_cpp.DMAT", L);
        }
    }
    //Debugging code
    // 
    // 
    //Eigen::MatrixXd err = Aeq * B_test;
    //std::cout << "error " << err <<  std::endl;
   // assert((Aeq * B_test).norm() < 1e-6);
}