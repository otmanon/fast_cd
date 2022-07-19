#include "QuasiNewtonSolver.h"
#include <igl/cat.h>
#include "line_search.h"
#include <iostream>
#include "augment_with_linear_constraints.h"
#include "igl/min_quad_with_fixed.h"
QuasiNewtonSolver::QuasiNewtonSolver()
{}


void QuasiNewtonSolver::precompute(const Eigen::SparseMatrix<double>& Q)
{
    Eigen::SparseMatrix<double> A;
    ldlt_precomp.compute(A);
    if (ldlt_precomp.info() != Success) {
        std::cout << "LLT factorization of system matrix failed" << std::endl;
        exit(0);
    }
}


void QuasiNewtonSolver::precompute_with_equality_constraints(const Eigen::SparseMatrix<double>& A,  const Eigen::SparseMatrix<double>& S)
{

    this->S = S;
    Eigen::SparseMatrix<double> Q;

    augment_with_linear_constraints(A, S, Q);
    // std::cout << Qeq.rowwise().sum() << std::endl;
    // ldlt_precomp.compute(Q);
    ldlt_precomp.compute(Q);
  

    
    if (ldlt_precomp.info() != Success) {
        std::cout << "LDLT factorization of system matrix failed" << std::endl;
        exit(0);
    }
}

Eigen::VectorXd QuasiNewtonSolver::solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, bool do_line_search, bool to_convergence, double max_iters)
{
    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g, ub;
    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows());
    rhs.setZero();

    z_next = z - S.transpose() * S * z;
    bool converged = false;
    double diff = 0;
    int i = 0;
    while (!converged)
    {
        alpha = 2.0;
        z_prev = z_next;
        g = grad_f(z_next);

       
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        dz = ldlt_precomp.solve(rhs);

        if (do_line_search)
        {
            const double energy0 = f(z_next);
            int line_search_step = 0;
            const double threshold = g.transpose() * dz.topRows(z_next.rows());
            do
            {
                alpha *= 0.5;
                z_next = z_prev + alpha * dz.topRows(z_next.rows());
                energy = f(z_next);
                //   printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e , threshold : %e\n", line_search_step, energy0, energy, alpha, threshold);
                line_search_step += 1;
            } while (energy > energy0 + 1e-5 && line_search_step < 100);
        }
        else
        {
            z_next += dz;
        }

        //stopping criteria
        i += 1;

        diff = (z_next - z_prev).norm();
     //   printf("iter %i , diff %g  \n ", i, diff);
        if (diff < 1e-5) //assuming unit height, can't really see motions on screen smaller than this value
        {
            converged = true;
        }
        else if (!to_convergence && i == max_iters)
        {
            break;
        }
        // std::cout << alpha << std::endl;
    }

    printf("Converged after %i iterations, with %g difference \n ", i, diff);

    return z_next;
}

Eigen::VectorXd QuasiNewtonSolver::solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0, bool do_line_search, bool to_convergence, double max_iters)
{
    assert(S.rows() == bc0.rows() && "Gave linear equality constraint rhs, but don't have a linear equality constraint matrix precomputed. \
                                                    Make sure you called precompute_with_constraints(Q, Qeq) before this");


    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g;

    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + S.rows()), dir;
    rhs.setZero();
    z_next = z;
    //project to constraint space first (this helps with convergence of high stiffness)s
    z_next = z + S.transpose() * (llt_proj.solve(bc0 - S * z));
   //  z_next = z - S.transpose() * (S * z - bc0);
    double e0 = 1;
    double e = e0 + 1;

    bool converged = false;
    double diff = 0;
    int i = 0;
    while(!converged)
    {
        e0 = e;
        z_prev = z_next;
        g = grad_f(z_next);

        //if (g.squaredNorm() < threshold)
        //     break;
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        rhs.bottomRows(S.rows()).setZero();
        const Eigen::VectorXd dir = ldlt_precomp.solve(rhs);
        dz = dir.topRows(z_next.rows());

        //z_next += dz; // set this for  now, though we should have a flag that skips line search
        if (do_line_search)
        {
            alpha = 2.0;
            double e0 = f(z_next);
            int line_search_step = 0;
            do
            {
                alpha *= 0.5;
                z_next = z_prev + alpha * dz;
                e = f(z_next);
              //  printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e \n", line_search_step, e0, e, alpha);
                line_search_step += 1;
            } while (e > e0 + 1e-9 && line_search_step < 100);
        }else
        {
            z_next += dz;
        }

        //stopping criteria
        i += 1;

        diff = (z_next - z_prev).norm();

        printf("iter %i , diff %g  \n ", i, diff);
        if (diff < 1e-5) //assuming unit height, can't really see motions on screen smaller than this value
        {
            converged = true;
        }
        else if (!to_convergence && i == max_iters)
        {
            break;
        }
    }
    printf("Converged after %i iterations, with %g difference \n ", i, diff);


    return z_next;
}