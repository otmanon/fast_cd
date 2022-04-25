#include "compute_modes.h"
#include <iostream>
#include <cassert>
#include <igl/sort.h>
#include <igl/get_seconds.h>
#include <igl/matlab/matlabinterface.h>

#ifdef FAST_CD_USE_SUITESPARSE
#include <Eigen/UmfPackSupport>
#include "UMFPACKSymShiftInvert.h"
#include "CholModSymShiftInvert.h"
#include <Eigen/CholmodSupport>
#include "CholModSimplicialLDLTSymShiftInvert.h"
#include "SimplicialLDLTSymShiftInvert.h"
#include <Spectra/SymGEigsShiftSolver.h>

using namespace Spectra;
#endif

void compute_modes(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, int r, Eigen::MatrixXd& U, Eigen::VectorXd& S)
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
  /*
    #ifdef FAST_CD_USE_SUITESPARSE
    using OpType = UMFPACKSymShiftInvert;
    using BopType = MatProd;

    BopType Bop(B);

    Eigen::UmfPackLU<Eigen::SparseMatrix<double>>* ldlt_solver = new  Eigen::UmfPackLU<Eigen::SparseMatrix<double>>();
   
    OpType op(A, B, ldlt_solver);
    SymGEigsShiftSolver<OpType, BopType, GEigsMode::ShiftInvert> geigs(op, Bop, r, r*2.0, 0);
  
    double total_start = igl::get_seconds();
    double t_start = igl::get_seconds();
    printf("Eigen::UMFPACKLU factorization beginning... \n");
  //   ldlt_solver->setMode(Eigen::CholmodMode::CholmodLDLt);
    op.set_mat(A, B);
    Bop.set_mat(B);
  //  op.set_solver(ldlt_solver);
    op.set_custom_shift(0);
    printf("Eigen::Eigen::UMFPACKLU succeeded in %g seconds ... \n", igl::get_seconds() - t_start);
 
  
    geigs.init();
    t_start = igl::get_seconds();
   // int nconv = geigs_umf.compute(SortRule::LargestMagn);
    std::cout << "Computing eigenvalues/eigenvectors using shift invert mode..." << std::endl;
   int nconv = geigs.compute(SortRule::LargestMagn);
    //B_spectra.resize(V.rows() * 3, 1);
    //B_spectra.setZero();
    if (geigs.info() == CompInfo::Successful)
    {
        printf("Found eigenvectors in %g seconds ... \n", igl::get_seconds() - t_start);
        printf("Total decomposition time %g... \n", igl::get_seconds() - total_start);
        U = geigs.eigenvectors();
      //  Eigen::VectorXd norms = (U.transpose() * B * U).diagonal();
      //  std::cout << norms << std::endl;
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
    assert((geigs.info() == CompInfo::Successful) && "Eigendecomposition not succeeded");
    #endif FAST_CD_USE_SUITESPARSE
    */
}
