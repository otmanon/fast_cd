#include "DenseQuasiNewtonSolver.h"
#include <igl/cat.h>
#include "line_search.h"
#include <iostream>
#include "augment_with_linear_constraints.h"
#include <igl/matlab/matlabinterface.h>
DenseQuasiNewtonSolver::DenseQuasiNewtonSolver()
{}


void DenseQuasiNewtonSolver::precompute(const Eigen::MatrixXd& Q)
{
    Eigen::MatrixXd A;
    llt_precomp.compute(Q);
    if (llt_precomp.info() != Success) {
        std::cout << "LLT factorization of system matrix failed" << std::endl;

        ldlt_precomp.compute(Q);
        if (ldlt_precomp.info() != Success)
        {
            std::cout << "ldlt failed too!" << std::endl;
        }
     //   exit(0);
    }
}


void DenseQuasiNewtonSolver::precompute_with_equality_constraints(const Eigen::MatrixXd& A, const Eigen::MatrixXd& S)
{

    this->S = S;
    Eigen::MatrixXd Q;


   // FullPivLU<MatrixXd> lu_decomp2(A);
   // std::cout << "rank of A : " << lu_decomp2.rank() << std::endl;
   // std::cout << "number of row of A : " << A.rows() << std::endl;
    augment_with_linear_constraints(A, S, Q);

  // FullPivLU<MatrixXd> lu_decomp1(S);
  // double constraint_rank_row = lu_decomp1.rank();
  // FullPivLU<MatrixXd> lu_decomp2(S.transpose());
  // double constraint_rank_col = lu_decomp2.rank();
  // FullPivLU<MatrixXd> lu_decomp3(S * S.transpose());
  // double constraint_rank = lu_decomp3.rank();
  // printf("row_rank %g  \n", constraint_rank_row);
  // printf("col_rank %g  \n", constraint_rank_col);
  // printf("rank %g  \n", constraint_rank);
  //  Engine* engine;
  //  Eigen::MatrixXd H = S * S.transpose();
  //  Eigen::MatrixXd I(S.rows(), S.rows());
  //  I.setIdentity();
  //  H += 1e-8 * I;

  // llt_proj.compute(H);
  // Eigen::VectorXd z, v;
  // Eigen::VectorXd ones = Eigen::VectorXd::Ones(H.rows());
  // if (H.rows() > 0)
  // {
  // z = llt_proj.solve(ones);
  // v = H * z;
  // std::cout << " Partial differentce error : " << (v - ones).norm() << std::endl;
  //  FullPivLU<MatrixXd> lu_decomp1(H);
  //  std::cout << "rank of H : " << lu_decomp1.rank() << std::endl;
  //  std::cout << "number of row of H : " << H.rows() << std::endl;
  //
  //  Eigen::MatrixXd testH = lu_decomp1.reconstructedMatrix();
  //  std::cout << "H difference " << (H - testH).norm() << std::endl;
  //
  //  Eigen::MatrixXd ST = S.transpose();
  //  FullPivLU<MatrixXd> lu_decomp4(ST);
  //  std::cout << "rank of S : " << lu_decomp4.rank() << std::endl;
  //  std::cout << "number of row of S : " << S.rows() << std::endl;

 //  }
 //  igl::matlab::mlinit(&engine);
 //  igl::matlab::mlsetmatrix(&engine, "S", S);
 //  igl::matlab::mleval(&engine, "spy(S)");
    // std::cout << Qeq.rowwise().sum() << std::endl;
    // ldlt_precomp.compute(Q);
   // Eigen::VectorXd rowsumH = H.rowwise().sum();
   // Eigen::VectorXd rowsumS = S.rowwise().sum();
    ldlt_precomp.compute(Q);
  //  ones = Eigen::VectorXd::Ones(Q.rows());
  //  z = lu_precomp.solve(ones);
  //  v = Q * z;
  //  Eigen::MatrixXd testQ = lu_precomp.reconstructedMatrix();
  //  std::cout << "Q difference " << (Q - testQ).norm() << std::endl;
  //  std::cout << " Full differentce error : " << (v - ones).norm() << std::endl;
  //
  //  std::cout << "Is invertible : " << lu_precomp.isInvertible() << std::endl;
  //
  //  std::cout << "rank of Q : " << lu_precomp.rank() << std::endl;
  //  std::cout << "number of row of Q : " << Q.rows() << std::endl;
    //llt_precomp.compute(A);
  // if (ldlt_precomp.info() != Success) {
  //     std::cout << "LDLT factorization of system matrix failed" << std::endl;
  //     exit(0);//
  // }
}

Eigen::VectorXd DenseQuasiNewtonSolver::solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, bool do_line_search, bool to_convergence, double  max_iters)
{
    Eigen::VectorXd z_prev, z_next;
    Eigen::VectorXd dz, g, ub;
    double alpha, energy, error = 1, threshold = 1e-9;

    Eigen::VectorXd rhs(z.rows());
    rhs.setZero();

    z_next = z;
    int i = 0;
    bool converged = false;
    double diff = 0;
    while (!converged)
    {
        alpha = 2.0;
        z_prev = z_next;
        g = grad_f(z_next);

      
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        dz = llt_precomp.solve(rhs);

        // z_next += dz;
          //itty bitty line search. should have a flag to do this instead of just stepping by 1. Not too costly tho, especially once we reduce
        if (do_line_search)
        {
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
            } while (energy > energy0 + 1e-5 && line_search_step < 100);
        }
        else
        {
            z_next += dz;
        }


        //stopping criteria
        i += 1;

        diff = (z_next - z_prev).norm();
        if (diff < 1e-6) //assuming unit height, can't really see motions on screen smaller than this value
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

Eigen::VectorXd DenseQuasiNewtonSolver::solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
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
    //project to constraint space first (this helps with convergence of high stiffness)
  //  z_next = -z + S.transpose() * (llt_proj.solve(bc0 + S * z));
    z_next = z + S.transpose() * (llt_proj.solve(bc0 - S * z));

    Eigen::VectorXd z_test = S * z;
    double e0 = 1;
    double e = e0 + 1;
    int i = 0;
    bool converged = false;
    double diff = 0;
    while (!converged)
    {
        e0 = e;
        z_next = z_next + S.transpose() * (llt_proj.solve(bc0 - S * z_next));
        z_test = S * z_next;
        std::cout << z_test.norm() << ": constraint enforcement" << std::endl;
        z_prev = z_next;
        g = grad_f(z_next);
        
     
        rhs.topRows(g.rows()) = -g;             //leave rhs dealing with constraints = 0 always
        rhs.bottomRows(S.rows()).setZero();
        const Eigen::VectorXd dir = ldlt_precomp.solve(rhs);
        dz = dir.topRows(z_next.rows());

       // Eigen::VectorXd v = lu_precomp.appl(rhs.transpose());
       // z_next += dz;

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
           //     printf("line_search_iter : %i, energy0 : %e, energy : %e,  alpha : %e \n", line_search_step, e0, e, alpha);
               line_search_step += 1;
           } while (e > e0 + 1e-9 && line_search_step <10);
        }
        else
        {
            z_next += dz;
        }

        
        //stopping criteria
        i += 1;
  
        diff = (z_next - z_prev).norm();
        if (diff < 1e-6) //assuming unit height, can't really see motions on screen smaller than this value
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