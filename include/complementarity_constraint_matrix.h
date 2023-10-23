#pragma once

#include "momentum_leaking_matrix.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/diag.h>

/*
Computes the complementarity constraint matrix for a mesh V, T with rig jacobian J, and optional input momentum leaking parameters
d

Input : 
V
T
J
(Optional)
D - n x 1 list of entries that fill the diagonal of the momentum leaking matrix
*/

using namespace Eigen;
void complementarity_constraint_matrix(const MatrixXd& V, const MatrixXi& T,
    const SparseMatrix<double>& J, const VectorXd& d, SparseMatrix<double>& Aeq)
{
    assert(d.rows() == V.rows() && "momentum leaking entries must have same dimensinoality \
                                    as vertices ");
    SparseMatrix<double> D(d.rows(), d.rows());
    igl::diag(d, D);
    D = igl::repdiag(D, V.cols());
    SparseMatrix<double>  M;
    igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    M = igl::repdiag(M, V.cols());
    Aeq = (J.transpose() * M * D);

}

void complementarity_constraint_matrix(const MatrixXd& V, const MatrixXi& T,
    const SparseMatrix<double>& J, SparseMatrix<double>& D, SparseMatrix<double>& Aeq)
{

    VectorXd d = D.diagonal();
    complementarity_constraint_matrix(V, T, J, d, Aeq);
}

void complementarity_constraint_matrix(const MatrixXd& V, const MatrixXi& T,
    const SparseMatrix<double>& J, std::string constraint_type, SparseMatrix<double>& Aeq)
{

    SparseMatrix<double>  D(V.rows(), V.rows());
    D.setIdentity();
    if (constraint_type == "momentum_leak")
    {
        momentum_leaking_matrix(V, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
    }
    complementarity_constraint_matrix(V, T, J, D, Aeq);
}