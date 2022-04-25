#include "CDNewtonSolver.h"
#include <igl/cat.h>
#include "line_search.h"
#include <iostream>
#include "augment_with_linear_constraints.h"
CDNewtonSolver::CDNewtonSolver(const int max_iter, const double tolerance, const const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq, double alpha, int max_iter_line_search) :
	max_iters(max_iter), tolerance(tolerance), alpha(alpha), max_iter_line_search(max_iter_line_search)
{
	Eigen::VectorXi ui;
	igl::min_quad_with_fixed_precompute(Q, ui, Qeq, false, precomp);
}

CDNewtonSolver::CDNewtonSolver(const int max_iter, const double tolerance, double alpha, int max_iter_line_search) :
    max_iters(max_iter), tolerance(tolerance), alpha(alpha), max_iter_line_search(max_iter_line_search)
{
    precomp;
}


void CDNewtonSolver::precompute(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq)
{
    Eigen::VectorXi ui;
   // igl::min_quad_with_fixed_precompute(Q, ui, Qeq, false, precomp);
    CompEq = Qeq;
    Eigen::SparseMatrix<double> A;
    augment_with_linear_constraints(Q, Qeq, A);
    ldlt_precomp.compute(A);
}


void CDNewtonSolver::precompute_with_constraints(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& Aeq, const Eigen::SparseMatrix<double>& Seq)
{
    Eigen::VectorXi bi;
    CompEq = Aeq;
    Eigen::SparseMatrix<double> Heq, Q;
    igl::cat(1, Aeq, Seq, Heq);
    Qeq = Heq;

    augment_with_linear_constraints(A, Heq, Q);

    //std::cout << Qeq.rowwise().sum() << std::endl;
    ldlt_precomp.compute(Q);

    if (ldlt_precomp.info() != Success) {
        std::cout << "LDLT factorization of system matrix failed" << std::endl;
        exit(0);
    }
 //   igl::min_quad_with_fixed_precompute(A, bi, Heq, false, precomp_constraints);

    //build Qeq matrix.  Should just be a Q.rows()x # clusters matrix, where each cluster column has a 1 at entry bI. Not too shabby guvna
   
    //std::vector<Eigen::Triplet<double>> tripletList;
    //tripletList.reserve(bi.rows());
    //for (int i = 0; i < bi.rows(); i++)
    //    tripletList.emplace_back(bi(i), i, 1);
    //
    //this->Qeq.resize(A.rows(), bi.rows());
    //this->Qeq.setFromTriplets(tripletList.begin(), tripletList.end());
    //
    //this->Qeq;
}

Eigen::VectorXd CDNewtonSolver::solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
     Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g, ub, Beq;
    Beq.resize(precomp.lagrange.rows()); Beq.setZero(); //complementary constraints
    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + CompEq.rows());
    rhs.setZero();
   
    z_next = z - CompEq.transpose() * CompEq * z;
    for (int i = 0; i < max_iters; i++)
    {
        alpha = 2.0;
        z_prev = z_next;
        g = grad_f(z_next);

        if (g.squaredNorm() < threshold)
            break;
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        dz = ldlt_precomp.solve(rhs);

      // z_next += dz;
        //itty bitty line search. should have a flag to do this instead of just stepping by 1. Not too costly tho, especially once we reduce
        const double energy0 = f(z_next);
       int line_search_step = 0;
       const double threshold =  g.transpose() * dz.topRows(z_next.rows());
       do
       {
           alpha *= 0.5;
           z_next = z_prev +  alpha*dz.topRows(z_next.rows());
           energy = f(z_next);
        //   printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e , threshold : %e\n", line_search_step, energy0, energy, alpha, threshold);
           line_search_step += 1;
       } while (energy > energy0 );
       // std::cout << alpha << std::endl;
    }
    
 
    return z_next;
}


Eigen::VectorXd CDNewtonSolver::solve_with_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0)
{
    assert(bc0.rows() <= Qeq.cols() && "Linear equality constraint matrix has different number of columns than constraints we are setting. Make sure we initialized it properly\
                                        by calling precompute_with_constraints");


    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g, Beq;
    Eigen::VectorXd lambda0 = Eigen::VectorXd::Zero(Qeq.rows());
    Eigen::VectorXd dlambda, lambda_prev, lambda_next=lambda0;
    
    Eigen::VectorXd r, r_new; //residuals used in line search
    
    Beq.resize(Qeq.rows()); Beq.setZero(); //complementary constraints
    Beq.bottomRows(bc0.rows()) = bc0;       //other constarints that may exist

    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + Qeq.rows()), dir ;
    rhs.setZero();
    z_next = z;
    //project to constraint space first (this helps with convergence of high stiffness)
    z_next = z_next - Qeq.transpose() * (Qeq * z_next - Beq);
    
    //
    double e0 = r.squaredNorm(), e;
    for (int i = 0; i < max_iters; i++)
    {
        alpha = 2.0;
     //   z_next = z_next - Qeq.transpose() * (Qeq * z_next - Beq);
        z_prev = z_next;
        lambda_prev = lambda_next;
        g = grad_f(z_next);
        
        Eigen::VectorXd test = Qeq.transpose() * lambda_next;
       // if (g.squaredNorm() < threshold)
       //     break;
        rhs.topRows(g.rows()) = -g - Qeq.transpose() * lambda_next;             //leave rhs dealing with constraints = 0 always
        rhs.bottomRows(Qeq.rows()) = -(Qeq  * z_next - Beq);
        const Eigen::VectorXd dir = ldlt_precomp.solve(rhs);
        dz = dir.topRows(z_next.rows());
        dlambda = dir.bottomRows(Qeq.rows());

        r = Eigen::VectorXd::Zero(rhs.rows());
        r.topRows(g.rows()) = grad_f(z_prev) + Qeq.transpose() * lambda_prev;
        r.bottomRows(lambda_prev.rows()) = (Qeq * z_prev- Beq);
       
        e0 = r.squaredNorm();
        int line_search_step = 0;
        do
        {
            alpha *= 0.5;
            z_next = z_prev + alpha * dz;
            lambda_next = lambda_prev +  alpha * dlambda;

            r_new = Eigen::VectorXd::Zero(rhs.rows());
            r_new.topRows(g.rows()) = grad_f(z_next) + Qeq.transpose() * lambda_next;
            r_new.bottomRows(lambda_prev.rows()) = (Qeq * z_next - Beq);
            e = r_new.squaredNorm();

           // printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e , threshold : %e\n", line_search_step, e0, e, alpha, threshold2);
            line_search_step += 1;
        } while (e > (1.0 + 1e-6)*e0  && line_search_step<10 );
        //std::cout << alpha << std::endl;
    }
    return z_next;
}