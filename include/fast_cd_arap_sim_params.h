#pragma once
#include "cd_sim_params.h"
#include "read_clusters_from_cache.h"
#include "read_modes_from_cache.h"
#include "lbs_jacobian.h"
#include "momentum_leaking_matrix.h"

#include <igl/massmatrix.h>

using namespace Eigen;
using namespace std;


struct fast_cd_arap_sim_params : cd_sim_params
{
		MatrixXd B;
		VectorXi labels;
		fast_cd_arap_sim_params() {};
		

		/*
		Contains all the parameters required to build a 
		fast Complementary Dynamics simulator.

		Inputs:
			X - n x 3 vertex geometry
			T - T x 4 tet indices
			B - 3n x m subspace matrix
			l - T x 1 clustering labels for each tet
			J - 3n x sol_p rig jacobian. If sim_constraint_type = none, this also 
									doubles as the equality constraint matrix
			mu - (double) first lame parameter
			lambda - (double) second lame parameter
			do_inertia - (bool) whether or not sim should have inertia 
					(if no, then adds Tik. regularizer to laplacian
			sim_constraint_type - (string) "none" or "cd" or "cd_momentum_leak"
		*/
		fast_cd_arap_sim_params(const MatrixXd& X, const MatrixXi& T, 
			const MatrixXd& B, const VectorXi&  l,
			const SparseMatrix<double>& J, SparseMatrix<double>& Aeq, double mu,
			double lambda, double h, bool do_inertia )
		{
			this->X = X;
			this->T = T;
			this->J = J;
			this->mu = mu * VectorXd::Ones(T.rows());
			this->lambda = VectorXd::Zero(T.rows());
			this->h = h;
			this->invh2 = 1.0 / (h * h);
			this->do_inertia = do_inertia;

			this->B = B;
			this->labels = l;
			this->Aeq = Aeq;
		}



		void set_equality_constraint(SparseMatrix<double>& Aeq)
		{
			this->Aeq = Aeq;
		}
	

};