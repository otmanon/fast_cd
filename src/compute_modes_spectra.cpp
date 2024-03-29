/*
#include <iostream>
#include <cassert>
#include <igl/sort.h>
#include <igl/get_seconds.h>
#include "compute_modes_spectra.h"

#include <Eigen/UmfPackSupport>
#include "UMFPACKSymShiftInvert.h"

#include <Spectra/SymGEigsShiftSolver.h>

using namespace Spectra;

void compute_modes_spectra(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, int r, Eigen::MatrixXd& U, Eigen::VectorXd& S)
{
        using OpType = UMFPACKSymShiftInvert;
        using BopType = MatProd;

        BopType Bop(B);

        Eigen::UmfPackLU<Eigen::SparseMatrix<double>>* ldlt_solver = new  Eigen::UmfPackLU<Eigen::SparseMatrix<double>>();

        OpType op(A, B, ldlt_solver);
        SymGEigsShiftSolver<OpType, BopType, GEigsMode::ShiftInvert> geigs(op, Bop, r, r*2.0, 0);

        double total_start = igl::get_seconds();
        double t_start = igl::get_seconds();

        printf("Eigen::UMFPACKLU factorization beginning... \n");
        op.set_mat(A, B);
        Bop.set_mat(B);
        op.set_custom_shift(0);
        printf("Eigen::UMFPACKLU factorization succeeded in %g seconds ... \n", igl::get_seconds() - t_start);

        geigs.init();
        t_start = igl::get_seconds();
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

        }else
        {
            printf("Mode computation failed!");
        }
        assert((geigs.info() == CompInfo::Successful) && "Eigendecomposition failed");
        
}
*/