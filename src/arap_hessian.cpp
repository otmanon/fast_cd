#include "arap_hessian.h"
#include <iostream>
#include <fstream>
#include "igl/repdiag.h"
#include "igl/cat.h"
#include "igl/cotmatrix.h"
#include "igl/arap_linear_block.h"
#include "igl/covariance_scatter_matrix.h"
using namespace igl;

void arap_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& H)
{
    const int dim = V.cols();
    int num_f = F.rows();
    int num_v = V.rows();
    //cotan laplacian and covariance scatter matrix
    Eigen::SparseMatrix<double> C, C_s, CSM;
    cotmatrix(V, F, C_s);
    C_s *= -3.0;
    
    repdiag(C_s, 3, C);
   


    covariance_scatter_matrix(V, F, ARAP_ENERGY_TYPE_ELEMENTS, CSM);

    //get covariacnes by doing CSM*repeat(V, dim, 1)
    Eigen::MatrixXd U = V.replicate(dim, 1);
    
    Eigen::MatrixXd S = CSM * U;                                          //stacked list of covariance matrices.
    S *= -3.0;
    Eigen::MatrixXd S_inv = Eigen::MatrixXd::Zero(S.rows(), S.cols());   //holds covariance matrices, with diagonal added on, inverted.
   

    Eigen::MatrixXd  cov(3, 3), cov_inv(3, 3);
    //build block diagonal covariance matrix
  
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(S.rows() * S.cols());
    for (int f = 0; f < F.rows(); f++)
    {
        // build this covariance matrix
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                cov(i, j) = S(i * num_f + f, j);
            }
        }

        Eigen::VectorXd diag = cov.diagonal();
        float diagsum = diag.sum();
        Eigen::MatrixXd d = diagsum*Eigen::MatrixXd::Identity(dim, dim);
        cov -=  d;
        cov_inv = cov.inverse();
        
        S_inv.block(dim*f, 0, dim , dim) = cov_inv;
        
        //there are now 9 triplet entries to fill out for this covariance matrix
        tripletList.emplace_back(T(num_f * 0 + f, num_f * 0 + f, cov_inv(0, 0)));
        tripletList.emplace_back(T(num_f * 0 + f, num_f * 1 + f, cov_inv(0, 1)));
        tripletList.emplace_back(T(num_f * 0 + f, num_f * 2 + f, cov_inv(0, 2)));

        tripletList.emplace_back(T(num_f * 1 + f, num_f * 0 + f, cov_inv(1, 0)));
        tripletList.emplace_back(T(num_f * 1 + f, num_f * 1 + f, cov_inv(1, 1)));
        tripletList.emplace_back(T(num_f * 1 + f, num_f * 2 + f, cov_inv(1, 2)));

        tripletList.emplace_back(T(num_f * 2 + f, num_f * 0 + f, cov_inv(2, 0)));
        tripletList.emplace_back(T(num_f * 2 + f, num_f * 1 + f, cov_inv(2, 1)));
        tripletList.emplace_back(T(num_f * 2 + f, num_f * 2 + f, cov_inv(2, 2)));
    }
    Eigen::SparseMatrix<double> SS(3 * num_f, 3 * num_f);
    SS.setFromTriplets(tripletList.begin(), tripletList.end());


    //Assemble B matrix
    Eigen::SparseMatrix<double> KX, KY, KZ, NKX, NKY, NKZ,Z(num_f, V.rows());
    arap_linear_block(V, F, 0, ARAP_ENERGY_TYPE_ELEMENTS, KX);
    arap_linear_block(V, F, 1, ARAP_ENERGY_TYPE_ELEMENTS, KY);
    arap_linear_block(V, F, 2, ARAP_ENERGY_TYPE_ELEMENTS, KZ);
    KX = KX.transpose().eval(); KY = KY.transpose().eval(); KZ = KZ.transpose().eval();
    NKX = -KX; NKY = -KY; NKZ = -KZ;
    Eigen::SparseMatrix<double> top_row = cat(2, cat(2, Z, NKZ), KY);
    Eigen::SparseMatrix<double> mid_row = cat(2, cat(2, KZ, Z), NKX);
    Eigen::SparseMatrix<double> bot_row = cat(2, NKY, cat(2, KX, Z));
    Eigen::SparseMatrix<double> B = -3.0*cat(1, cat(1, top_row, mid_row), bot_row);


    H = C - B.transpose() * SS.transpose() * B;

   // Eigen::SparseMatrix<double> B = SS.transpose();
   // float norm = (A - B).norm();

   // auto norm = (A - SS).norm();
    //assert((SS - SS.transpose()).norm() < 1e-8); //makes as
    //const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

    //std::ofstream file("S.csv");
    //if (file.is_open())
    //{
    //    file << S.format(CSVFormat) << '\n';
    //    file.close();
    //}

    //std::ofstream file2("S_inv.csv");
    //if (file2.is_open())
    //{
    //    file2 << S_inv.format(CSVFormat) << '\n';
    //    file2.close();
    //}
   // std::cout << R << std::endl;


}