#include "compute_compressed_vibration_modes.h"

#include <igl/active_set.h>
#include <igl/cat.h>
#include <igl/repdiag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/matlab/matlabinterface.h>

//#include <igl/mosek/mosek_quadprog.h>

void compute_compressed_vibration_modes(Eigen::SparseMatrix<double>& H, Eigen::SparseMatrix<double>& M, int num_modes, Eigen::MatrixXd& U, double mu, bool to_convergence , int max_iters){
	Eigen::SparseMatrix<double> Z;
	Z.resize(H.rows(), H.rows());
	Z.setZero();
	Eigen::SparseMatrix<double> Hn, Mn;
	Eigen::SparseMatrix<double> I;
	I.resize(H.rows(), H.rows());
	I.setIdentity();
	//H += 1e-8 * I;
	//H *= -1.0;
	Hn = -H; Mn = -M;


	// Matlab instance
	Engine* engine;
	igl::matlab::mlinit(&engine);
	


	Eigen::SparseMatrix<double> H_exp, M_exp, M_flat;
	
	Eigen::VectorXd M1 = M * Eigen::VectorXd::Ones(M.rows()) ;

	H_exp = igl::cat(1, igl::cat(2, H, Hn),
					igl::cat(2, Hn, H));

	M_exp = igl::cat(1, igl::cat(2, M, Z),
		igl::cat(2, Z, M));

	Eigen::VectorXd m = M_exp.diagonal();

	Eigen::VectorXd M1M1 = mu* igl::cat(1, M1, M1);

	Eigen::VectorXi known; Eigen::VectorXd bc;
	U.conservativeResize(H.rows(), 0);

	Eigen::MatrixXd U_exp;
	Eigen::VectorXd Beq, Bineq;
	Eigen::SparseMatrix<double> Aeq,  Aineq;

	Aineq = -igl::cat(1, igl::cat(2, I, Z),
						igl::cat(2, Z, I));
	Bineq.resize(Aineq.rows());
	Bineq.setZero();

	igl::matlab::mlsetmatrix(&engine, "H",  H);

	igl::matlab::mlsetmatrix(&engine, "M", M);
	igl::matlab::mlsetmatrix(&engine, "A", H_exp);
	igl::matlab::mlsetmatrix(&engine, "B", (Eigen::MatrixXd) M1M1);

	igl::matlab::mlsetmatrix(&engine, "Aieq", Aineq);
	igl::matlab::mlsetmatrix(&engine, "Bieq", (Eigen::MatrixXd) Bineq);
	// igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
	//igl::matlab::mlgetmatrix(&engine, "EV", U);
	//Eigen::VectorXd ux = Eigen::VectorXd::Zero(H_exp.rows());

	for (int i = 0; i < num_modes; i++)
	{
		int iter = 0;
		double diff = 1;
		Eigen::VectorXd ci = Eigen::VectorXd::Random(H.rows());
		ci.array() = ci.array() / ci.norm(); //normalized initial guess

		U.conservativeResize(U.rows(), U.cols() + 1);
		U.col(U.cols() - 1) = ci;
		U_exp.resize(U.rows() * 2, U.cols());
		U_exp.topRows(U.rows()) = U;
		U_exp.bottomRows(U.rows()) = -U;

		Beq.conservativeResize(Beq.rows() + 1);
		Beq.setZero();
		Beq(Beq.rows() - 1) = 1;


		Aeq = (M_exp * U_exp).sparseView().transpose();
		igl::matlab::mlsetmatrix(&engine, "Aeq", Aeq);
		igl::matlab::mlsetmatrix(&engine, "Beq", (Eigen::MatrixXd)Beq);
		while (true)
		{
			igl::active_set_params params;
			params.Auu_pd = false;
			Eigen::MatrixXd z;
		//	igl::SolverStatus status = igl::active_set(H_exp, M1M1,  known, bc, Aeq, Beq, Aineq, Bineq, lx, ux, params, z);
			
			igl::matlab::mleval(&engine, "ut = quadprog(A, B', Aieq, Bieq', Aeq, Beq)");
			Eigen::MatrixXd ut;
			igl::matlab::mlgetmatrix(&engine, "ut", z);

			Eigen::MatrixXd test;
		
			igl::min_quad_with_fixed_data<double> data;
			//Aeq.resize(0, 0);
			
			//igl::min_quad_with_fixed_precompute(H_exp, known, Aeq, false, data);
		//	Eigen::VectorXd z_test;
		//	igl::min_quad_with_fixed(H_exp, M1M1, known, bc, Aeq, Beq, false, z_test);
			Eigen::VectorXd u = z.topRows(H.rows()) - z.bottomRows(H.rows());
			ci = u.array() / u.norm();

			U.col(U.cols() - 1) = ci;
			U_exp.resize(U.rows() * 2, U.cols());
			U_exp.topRows(U.rows()) = U;
			U_exp.bottomRows(U.rows()) = -U;
			Aeq = (M_exp * U_exp).sparseView().transpose();


			

			igl::matlab::mlsetmatrix(&engine, "Aeq", Aeq);
			//igl::matlab::mlsetmatrix(&engine, "Beq", (Eigen::MatrixXd)Beq);
			iter += 1;
			//if (to_convergence) //not implemented
			//{
			//	if (diff < 1e-6)
			//		break;
			//}
			//else
			{
				if (iter == max_iters)
					break;
			}
		}
	}
}