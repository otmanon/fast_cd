#include "FastCDNewtonSolver.h"
#include "augment_with_linear_constraints.h"
#include "line_search.h"
FastCDNewtonSolver::FastCDNewtonSolver(const int max_iter, const double tolerance, const Eigen::MatrixXd& Q, double alpha, int max_iter_line_search) :
    max_iters(max_iter), tolerance(tolerance), alpha(alpha), max_iter_line_search(max_iter_line_search)
{
    Eigen::VectorXi ui;
    precomp = Eigen::LLT<Eigen::MatrixXd>();
    ldlt_precomp = Eigen::LDLT<Eigen::MatrixXd>();
    precomp.compute(Q);
}

FastCDNewtonSolver::FastCDNewtonSolver(const int max_iter, const double tolerance, double alpha, int max_iter_line_search) :
    max_iters(max_iter), tolerance(tolerance), alpha(alpha), max_iter_line_search(max_iter_line_search)
{
    precomp = Eigen::LLT<Eigen::MatrixXd>();
    ldlt_precomp = Eigen::LDLT<Eigen::MatrixXd>();
}

void FastCDNewtonSolver::precompute(const Eigen::MatrixXd& Q)
{
    Eigen::VectorXi ui;
    Eigen::SparseMatrix<double> Qeq;
    precomp.compute(Q);

    if (precomp.info() != Success) {
        std::cout << "LLT factorization of system matrix failed... maybe not psd?" << std::endl;
        exit(0);
    }
}

void FastCDNewtonSolver::precompute_with_constraints(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& Qeq)
{
    this->Qeq = Qeq;
    augment_with_linear_constraints(Q, Qeq, A);

    ldlt_precomp.compute(A);

    if (ldlt_precomp.info() != Success) {
        std::cout << "LDLT factorization of system matrix failed" << std::endl;
        exit(0);
    }

}

Eigen::VectorXd FastCDNewtonSolver::solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g;
    double alpha, e, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows());
    rhs.setZero();

    z_next = z;
    threshold = 1e-8;
    for (int i = 0; i < max_iters; i++)
    {
        alpha = 2.0;
        z_prev = z_next;
        g = grad_f(z_next);

        if (g.squaredNorm() < threshold)
            break;
        rhs = -g;             //leave rhs dealing with constraints = 0 always
        dz = precomp.solve(rhs);

      //  z_next += dz;
            //itty bitty line search. should have a flag to do this instead of just stepping by 1. Not too costly tho, especially once we reduce
        const double e0 = rhs.squaredNorm();
        int line_search_step = 0;
        const double threshold = abs(g.transpose() * dz);
        do
        {
            alpha *= 0.5;
            z_next = z_prev + alpha * dz;
            Eigen::VectorXd r = grad_f(z_next);
            e = r.squaredNorm();
            line_search_step += 1;
        } while (e > e0 );
      //  std::cout << alpha << std::endl;
    }
    return z_next;
}


Eigen::VectorXd FastCDNewtonSolver::solve_with_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f,  const Eigen::VectorXd& bc0)
{
    assert(Qeq.rows() == bc0.rows() && "Gave linear equality constraint rhs, but don't have a linear equality constraint matrix precomputed. \
                                                    Make sure you called precompute_with_constraints(Q, Qeq) before this");


    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g;
    Eigen::VectorXd lambda0 = Eigen::VectorXd::Zero(Qeq.rows());
    Eigen::VectorXd dlambda, lambda_prev, lambda_next = lambda0;

    Eigen::VectorXd r, r_new; //residuals used in line search

  
    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows() + Qeq.rows()), dir;
    rhs.setZero();
    z_next = z;
    //project to constraint space first (this helps with convergence of high stiffness)
    z_next = z_next - Qeq.transpose() * (Qeq * z_next - bc0);

    double e0 = r.squaredNorm(), e;
    for (int i = 0; i < max_iters; i++)
    {
        alpha = 2.0;
        z_prev = z_next;
        lambda_prev = lambda_next;
        g = grad_f(z_next);

        Eigen::VectorXd test = Qeq.transpose() * lambda_next;
        if (g.squaredNorm() < threshold)
             break;
        rhs.topRows(g.rows()) = -g - Qeq.transpose() * lambda_next;             //leave rhs dealing with constraints = 0 always
        rhs.bottomRows(Qeq.rows()) = -(Qeq * z_next - bc0);
        const Eigen::VectorXd dir = ldlt_precomp.solve(rhs);
        dz = dir.topRows(z_next.rows());
        dlambda = dir.bottomRows(Qeq.rows());

        r = Eigen::VectorXd::Zero(rhs.rows());
        r.topRows(g.rows()) = grad_f(z_prev) + Qeq.transpose() * lambda_prev; 
        r.bottomRows(lambda_prev.rows()) = (Qeq * z_prev - bc0);

        e0 = r.squaredNorm();
        int line_search_step = 0;
        do
        {
            alpha *= 0.5;
            z_next = z_prev + alpha * dz;
            lambda_next = lambda_prev + alpha * dlambda;

            r_new = Eigen::VectorXd::Zero(rhs.rows());
            r_new.topRows(g.rows()) = grad_f(z_next) + Qeq.transpose() * lambda_next;
            r_new.bottomRows(lambda_prev.rows()) = (Qeq * z_next - bc0);
            e = r_new.squaredNorm();

            // printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e , threshold : %e\n", line_search_step, e0, e, alpha, threshold2);
            line_search_step += 1;
        } while (e > (1.0 + 1e-6) * e0);
        //std::cout << alpha << std::endl;
    }
    return z_next;
}