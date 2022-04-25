#include "EigenHook.h"
#include "arap_hessian.h"
#include "igl/repdiag.h"
#include "igl/cat.h"
#include "compute_modes.h"
#include "augment_with_linear_constraints.h"
#include "momentum_leaking_matrix.h"
#include <igl/per_face_normals.h>
#include "igl/readMSH.h"
#include <igl/boundary_facets.h>
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>


void EigenHook::initSimulation()
{
    Eigen::VectorXi tritag, tettag;
    std::string data_dir = "../data/elephant/";
    std::string cache_dir = data_dir + "cache/";
    std::string B_file_dir = cache_dir + "B_affine_rig.eigcache";
    igl::readMSH("../data/elephant/elephant.msh", origQ, F, T, tritag, tettag);

    origQ = origQ.rowwise() - origQ.colwise().mean();
    origQ /= origQ.col(1).maxCoeff();

    igl::boundary_facets(T, F);
    Eigen::MatrixXi Ftempy = F.col(0);
    F.col(0) = F.col(1);

    F.col(1) = Ftempy;
    igl::per_face_normals(origQ, F, N_faces);
    Q = origQ;
    V.resizeLike(origQ);
    V.setZero();
    rig = AffineRig(origQ);
    
    Eigen::SparseMatrix<double> C, M, H;
    arap_hessian(Q, T, H);
    igl::cotmatrix(Q, T, C);
    C = igl::repdiag(C, 3);
    Eigen::SparseMatrix<double> I(H.rows(), H.cols());
    I.setIdentity();
    H = -H; // +1e-6 * I;
    C = C + 1e-6 * I;
    igl::massmatrix(Q, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    M = igl::repdiag(M, 3);
    Eigen::VectorXd m = M.diagonal();
    Eigen::VectorXd delta =1e-12*Eigen::VectorXd::Ones(rig.get_J().cols());
    Eigen::VectorXd diag = igl::cat(1, m, delta);
    int n = origQ.rows() * origQ.cols() + rig.get_J().cols();
    Eigen::MatrixXd M_exp = (diag.asDiagonal());
    Eigen::SparseMatrix<double> M_e = M_exp.sparseView();
   // M_exp.diagonal() = diag;
    Eigen::SparseMatrix<double> H_exp;
    Eigen::SparseMatrix<double> D;
    momentum_leaking_matrix(origQ, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
    D = igl::repdiag(D, 3);


    
    Eigen::SparseMatrix<double> A_eq = (rig.get_J().transpose() *(M *D));
    augment_with_linear_constraints(H, A_eq, H_exp);

    //  igl::eigs(H, M, 10, igl::EIGS_TYPE_SM, B_igl, S);
    Eigen::VectorXd S;//eigenvalues
    if (! read_binary(B_file_dir, B_spectra))
    {
        compute_modes(H_exp, M_e, r, B_spectra, S);
        write_binary(B_file_dir, B_spectra);
    }
    else
    {
        std::cout << "Using Cached Modes..." << std::endl;
    }

    // write_binary("Q_test.txt", Q);

    mode = 0;
    step = 0;
    dt = 1e-4;
    k = 1e-5;
}

bool EigenHook::simulateOneStep()
{
    bool use_spectra = true;
    //Q += dt * Eigen::MatrixXd::Ones(Q.rows(), Q.cols());
    step += 1;
    if (step % period == 0)
    {
        mode += 1;
        mode = mode % B_spectra.cols();
    }
    if (step - period * (step / period) == half_period)
    {
        scale *= -1; 
    }
      
    Eigen::MatrixXd B_flat = B_spectra.col(mode);
    Eigen::Map<Eigen::MatrixXd> B_mat(B_flat.data(), V.rows(), 3);
    Eigen::MatrixXd disp = scale * dt * B_mat;

    Q = Q + disp;

    return false;
}