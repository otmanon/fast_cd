#pragma once
#include "cd_sim_params.h"
#include "read_clusters_from_cache.h"
#include "read_modes_from_cache.h"
#include "lbs_jacobian.h"
#include "momentum_leaking_matrix.h"

#include <igl/massmatrix.h>

using namespace Eigen;
using namespace std;


struct fast_cd_arap_sim_params 
{

	//geometry n x dim
	MatrixXd X;

	//tet indices t x 3 (4)
	MatrixXi T;

	//timestep (s)
	double h;
	double invh2;

	// first lame parameter per tet (for ARAP it's basically young's modulus /2)
	VectorXd mu;

	//Equality constraints given as input from user
	SparseMatrix<double> Aeq;

	// linear rig jacobian! n*dim x b*dim*(dim+1)
	SparseMatrix<double> J;

	//whether to activate inertia or not. if false, will add a bit of regularization to quadratic term
	bool do_inertia;

	//Extra user defined quadratic term (like for mass springs or penalty constraints)
	SparseMatrix<double> Q;

	//subspace matrix for displacement (dn x m)
	MatrixXd B;

	// cluster labels for each tet  (l x 1)
	VectorXi labels;

		

	/*
	Contains all the parameters required to build a 
	fast Complementary Dynamics simulator.

	Inputs:
		X - n x d vertex geometry
		T - t x 3 (4) tet indices
		B - dn x m subspace matrix
		l - t x 1 clustering labels for each tet
		J - dn x p rig jacobian. If sim_constraint_type = none, this also 
								doubles as the equality constraint matrix

	    Aeq - c x dn linear equality constraint matrix
		mu - (double) first lame parameter
		lambda - (double) second lame parameter
		do_inertia - (bool) whether or not sim should have inertia 
				(if no, then adds Tik. regularizer to laplacian
	*/
	fast_cd_arap_sim_params(const MatrixXd& X, const MatrixXi& T, 
		const MatrixXd& B, const VectorXi&  l,
		const SparseMatrix<double>& J, const SparseMatrix<double>& Aeq, double mu,
		 double h, bool do_inertia )
	{
		this->X = X;
		this->T = T;
		this->J = J;
		this->mu = mu * VectorXd::Ones(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;

		this->B = B;
		this->labels = l;
		this->Aeq = Aeq;
	}


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
		Aeq - c x dn linear equality constraint matrix
		mu - T x 1 first lame parameter
		do_inertia - (bool) whether or not sim should have inertia
				(if no, then adds Tik. regularizer to laplacian
	*/
	fast_cd_arap_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		const SparseMatrix<double>& J, const SparseMatrix<double>& Aeq, 
		const VectorXd& mu, double h, bool do_inertia)
	{
		this->X = X;
		this->T = T;
		this->J = J;
		this->mu = mu;
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