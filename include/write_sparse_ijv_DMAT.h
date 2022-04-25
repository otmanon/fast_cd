#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/find.h>
#include <igl/writeDMAT.h>
#include <igl/cat.h>

bool write_sparse_ijv_DMAT(std::string filename, Eigen::SparseMatrix<double> A)
{
    Eigen::VectorXi I, J;
    Eigen::MatrixXi Ii, Ji;
    Eigen::VectorXd V;
    Eigen::MatrixXd Id, Jd, Vd;
    igl::find(A, I, J, V);
    Id = Eigen::Map<Eigen::MatrixXi>(I.data(), I.rows(), 1).cast<double>();
    Jd = Eigen::Map<Eigen::MatrixXi>(J.data(), J.rows(), 1).cast<double>();
    Vd = Eigen::Map<Eigen::MatrixXd>(V.data(), V.rows(), 1);
    Eigen::MatrixXd ijv = igl::cat(2, Id, (igl::cat(2, Jd, Vd)));
    return igl::writeDMAT(filename, ijv);
}