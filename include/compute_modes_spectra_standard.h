#include <iostream>
#include <cassert>
#include <igl/sort.h>
#include <igl/get_seconds.h>
#include "compute_modes_spectra.h"

#include <Spectra/SymGEigsShiftSolver.h>



using namespace Spectra;


void compute_modes_spectra_standard(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, int r, Eigen::MatrixXd& U, Eigen::VectorXd& S)
{

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;
    OpType op(A, B);
    BOpType Bop(B);

    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert>
        geigs(op, Bop, r, 2*r, 0.0);
    double total_start = igl::get_seconds();
    double t_start = igl::get_seconds();
    printf("Eigen::SparseLU factorization beginning... \n");
  
  

    geigs.init();
    printf("Eigen::SparseLU factorization succeeded in %g seconds ... \n", igl::get_seconds() - t_start);

    std::cout << "Computing eigenvalues/eigenvectors using shift invert mode..." << std::endl;
    int nconv = geigs.compute(SortRule::LargestMagn);

    if (geigs.info() == CompInfo::Successful)
    {
        printf("Found Spectra eigenvectors in %g seconds ... \n", igl::get_seconds() - t_start);
        printf("Total Spectra decomposition time %g... \n", igl::get_seconds() - total_start);
        U = geigs.eigenvectors();

        S = geigs.eigenvalues();
        Eigen::MatrixXd S_mat = Eigen::MatrixXd(S.cwiseAbs());
        //sort these according to S.abs()
        Eigen::MatrixXi I;
        Eigen::MatrixXd S_sorted;
        igl::sort(S_mat, 1, true, S_sorted, I);

        Eigen::MatrixXd U_sorted = Eigen::MatrixXd::Zero(U.rows(), U.cols());
        for (int i = 0; i < S.rows(); i++)
        {
            U_sorted.col(i) = U.col(I(i));
        }

        S = S_sorted;
        U = U_sorted;

    }
    else
    {
        printf("Mode computation failed!");
    }
    assert((geigs.info() == CompInfo::Successful) && "Eigendecomposition not succeeded");

}
