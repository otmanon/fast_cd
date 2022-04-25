#include "CDLocalGlobalSolver.h"
#include <iostream>
CDLocalGlobalSolver::CDLocalGlobalSolver(
    int max_iters, double tol, const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq)
    : max_iters(max_iters), tol(tol)
{
    Eigen::VectorXi ui;
    igl::min_quad_with_fixed_precompute(Q, ui, Qeq, false, precomp);
};

CDLocalGlobalSolver::CDLocalGlobalSolver(
    int max_iters, double tol)
    : max_iters(max_iters), tol(tol)
{
    precomp;
};

void CDLocalGlobalSolver::precompute(
    const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq)
{
     Eigen::VectorXi ui;
     igl::min_quad_with_fixed_precompute(Q, ui, Qeq, false, precomp);
};

Eigen::VectorXd CDLocalGlobalSolver::solve(const Eigen::VectorXd& z,
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& local_step,
    std::function<Eigen::VectorXd(const Eigen::MatrixXd&, const igl::min_quad_with_fixed_data<double>&)>&global_step)
{
    double res = 1;
    Eigen::VectorXd z0, z_next;
    Eigen::MatrixXd R_stack;
    z_next = z;
    int i = 0;
    for (i = 0; i < max_iters; i++)
    {
        z0 = z_next;
        //ger best fit rotations with current reduced coefficients z
        R_stack = local_step(z_next);
        //global step will return displacement but will be split up into different cases if we are using complementarity or not
        z_next = global_step(R_stack, precomp);
        res = (z_next - z0).squaredNorm();

        if (res < tol)
            break;
    }
    return z_next;
}



