#include "DenseQuasiNewtonSolver.h"
#include <igl/cat.h>
#include "line_search.h"
#include <iostream>
#include "augment_with_linear_constraints.h"

DenseQuasiNewtonSolver::DenseQuasiNewtonSolver(const int max_iter, const double tolerance, double alpha, int max_iter_line_search) :
    max_iters(max_iter), tolerance(tolerance), alpha(alpha), max_iter_line_search(max_iter_line_search)
{}


void DenseQuasiNewtonSolver::precompute(const Eigen::MatrixXd& Q)
{
    Eigen::MatrixXd A;
    ldlt_precomp.compute(Q);
    if (ldlt_precomp.info() != Success) {
        std::cout << "LLT factorization of system matrix failed" << std::endl;
        exit(0);
    }
}


void DenseQuasiNewtonSolver::precompute_with_equality_constraints(const Eigen::MatrixXd& A, const Eigen::MatrixXd& S)
{

    this->S = S;
    Eigen::MatrixXd Q;

    augment_with_linear_constraints(A, S, Q);

    llt_proj.compute(S * S.transpose());
    // std::cout << Qeq.rowwise().sum() << std::endl;
    // ldlt_precomp.compute(Q);
    ldlt_precomp.compute(Q);
    //llt_precomp.compute(A);
   if (ldlt_precomp.info() != Success) {
       std::cout << "LDLT factorization of system matrix failed" << std::endl;
       exit(0);
   }
}

Eigen::VectorXd DenseQuasiNewtonSolver::solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g, ub;
    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows());
    rhs.setZero();

    z_next = z - S.transpose() * S * z;

    for (int i = 0; i < max_iters; i++)
    {
        alpha = 2.0;
        z_prev = z_next;
        g = grad_f(z_next);

        if (g.squaredNorm() < 1e-5)
            break;
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        dz = ldlt_precomp.solve(rhs);

        // z_next += dz;
          //itty bitty line search. should have a flag to do this instead of just stepping by 1. Not too costly tho, especially once we reduce
        const double energy0 = rhs.norm(); // f(z_next);
        int line_search_step = 0;
        const double threshold = g.transpose() * dz.topRows(z_next.rows());
        do
        {
            alpha *= 0.5;
            z_next = z_prev + alpha * dz.topRows(z_next.rows());
            energy = grad_f(z_next).norm();
            //   printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e , threshold : %e\n", line_search_step, energy0, energy, alpha, threshold);
            line_search_step += 1;
        } while (energy > energy0 + 1e-5 && line_search_step < 10);
        // std::cout << alpha << std::endl;
    }
    std::cout << "energy : " << f(z_next) << std::endl;
    return z_next;
}

Eigen::VectorXd DenseQuasiNewtonSolver::solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0)
{
   
    assert(S.rows() == bc0.rows() && "Gave linear equality constraint rhs, but don't have a linear equality constraint matrix precomputed. \
                                                    Make sure you called precompute_with_constraints(Q, Qeq) before this");


    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g;

    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + S.rows()), dir;
    rhs.setZero();
    z_next = z;
    //project to constraint space first (this helps with convergence of high stiffness)
  //  z_next = -z + S.transpose() * (llt_proj.solve(bc0 + S * z));
    z_next = z + S.transpose() * (llt_proj.solve(bc0 - S * z));
    double e0 = 1;
    double e = e0 + 1;
    for (int i = 0; i < max_iters; i++)
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
       // z_next += dz;
        
       alpha = 2.0;
       double e0 = f(z_next);
       int line_search_step = 0;
       do
       {
           alpha *= 0.5;
           z_next = z_prev + alpha * dz;          
           e = f(z_next);
       //     printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e \n", line_search_step, e0, e, alpha);
           line_search_step += 1;
       } while (e > e0 + 1e-9 && line_search_step <10);
    }
    return z_next;
    
    /*
    assert(bc0.rows() == S.rows() && "Linear equality constraint matrix has different number of rows than rhs constraints we are setting. Make sure we initialized it properly\
                                        by calling precompute_with_equality_constraints");
    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, dlambda, g, Beq;

    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + S.rows()), dir,r , r_new;
    rhs.setZero();
    z_next = z;
    //project back to constraint space first (this helps with convergence of high stiffness)
    z_next = z_next - S.transpose() * (S * z_next -  bc0);

    Eigen::VectorXd test = S * z_next;
    double c, m, e0, e;
    c = 0.3;
    int line_search_step = 0;
    Eigen::VectorXd lambda_prev = Eigen::VectorXd::Zero(S.rows());
    Eigen::VectorXd lambda_next = lambda_prev;
    for (int i = 0; i < max_iters; i++)
    {
        z_prev = z_next;
        lambda_prev = lambda_next;
        g = grad_f(z_prev);

      //  if (g.squaredNorm() < 1e-5)
      //     break;
        rhs.topRows(g.rows()) = -g -S.transpose() * lambda_prev;              //leave rhs dealing with constraints = 0 always
        rhs.bottomRows(S.rows()) = bc0 - S * z_prev;
        const Eigen::VectorXd dir = ldlt_precomp.solve(rhs);
        dz = dir.topRows(z_prev.rows());
        dlambda = dir.bottomRows(S.rows());
        //z_next += dz;
      //  lambda += dlambda;
  
        //itty bitty line search

        r = Eigen::VectorXd::Zero(rhs.rows());
        r.topRows(g.rows()) = grad_f(z_prev) + S.transpose() * lambda_prev;
        r.bottomRows(lambda_prev.rows()) = (S* z_prev - bc0);

        alpha = 2.0;
        line_search_step = 0;
        e0 = rhs.squaredNorm();// f(z_next);
        do
        {
            alpha *= 0.5;
            z_next = z_prev + alpha * dz;
            lambda_next = lambda_prev + alpha * dlambda;
            g = grad_f(z_next);
            r_new = Eigen::VectorXd::Zero(rhs.rows());
            r_new.topRows(g.rows()) = grad_f(z_next) + S.transpose() * lambda_next;
            r_new.bottomRows(lambda_prev.rows()) = (S * z_next - bc0);
            e = r_new.squaredNorm();
            printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e \n", line_search_step, e0, e, alpha);
            line_search_step += 1;
        } while (e > (1.0 + 1e-5)*e0 );
        
        //printf("line_search_iters : %i, energy : %e,  alpha : %e,  grad_norm : %e, dir_norm : %e\n", line_search_step, e, alpha, dz.norm(), dz.norm());
        //std::cout << alpha << std::endl;
    }
    //  std::cout << "energy : " << f(z_next) << std::endl;
    return z_next;
    */
    
}