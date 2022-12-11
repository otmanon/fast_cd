#pragma once
#include "fast_cd_arap_sim.h"
#include "lbs_jacobian.h"
#include "selection_matrix.h"
#include "sim_state.h"
#include "fast_ik_arap_local_global_solver.h"
#include <igl/cat.h>
struct fast_ik_sim : fast_cd_arap_sim
{

	fast_ik_sim() {};
	/*
	Builds a FAST IK [Jacobson et. al 2014] Simulator,
	as a wrapper of FAST CD
	Inputs :
		V - (n x 3) Mesh Vertices
		T - (T x 4) tet indices
		W - n x m rig weights
		l - n x 1 cluster labels
		bI - bI x 1 constraint indices into V
		cd_arap_local_global_solver_params - (nuff said)
	*/
	fast_ik_sim(const MatrixXd& V, const MatrixXi& T, const MatrixXd& W,
		const VectorXi& l, const VectorXi& bI,
		const cd_arap_local_global_solver_params& solver_params)
	{
		SparseMatrix<double> B;
		lbs_jacobian(V, W, B);

		MatrixXd Bd = B.toDense();

		SparseMatrix<double> S;
		selection_matrix(bI, V.rows(), V.cols(), S);

		SparseMatrix<double> J(0, B.rows());

		params = new fast_cd_sim_params(V, T,
			Bd, l, J, S, 1, 0, 1e-2, false);

		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp(*((fast_cd_sim_params*)params));

		sp = fcd_sp; // fast_cd_arap_static_precomp(sim_params);
		dp = new fast_cd_arap_dynamic_precomp();

		sol = new fast_ik_arap_local_global_solver(fcd_sp->BAB, fcd_sp->AeqB, solver_params);

	};


	/*
	Builds a FAST IK [Jacobson et. al 2014] Simulator,
	as a wrapper of FAST CD
	Inputs :
		V - (n x 3) Mesh Vertices
		T - (T x 4) tet indices
		W - n x m rig weights
		l - n x 1 cluster labels
		S - c x dim n equality constraint matrix
		cd_arap_local_global_solver_params - (nuff said)
	*/
	fast_ik_sim(const MatrixXd& V, const MatrixXi& T, const MatrixXd& W,
		const VectorXi& l, const SparseMatrix<double>& S,
		const cd_arap_local_global_solver_params& solver_params)
	{
		SparseMatrix<double> B;
		lbs_jacobian(V, W, B);

		MatrixXd Bd = B.toDense();

		//no rig in these dynamics
		SparseMatrix<double> J(B.rows(), 0);


		params = new fast_cd_sim_params(V, T,
			Bd, l, J, S, 1, 0, 1e-2, false);

		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp(*((fast_cd_sim_params*)params));

		sp = fcd_sp; // fast_cd_arap_static_precomp(sim_params);
		dp = new fast_cd_arap_dynamic_precomp();

		sol = new fast_ik_arap_local_global_solver(fcd_sp->BAB, fcd_sp->AeqB, solver_params);
	};

	/*
	Steps the FAST IK simulation
	Inputs :
	z - (m x 1) initial guess for z
	state - simulation state (importatn for inertia if turned on)
	f_ext - (m x 1) ecternal force
	bc - (c x 1) DISPLACEMENT constraints on the simulation (must mattch with params.Aeq.rows())
	Outputs:
	z_next - (m x 1) next timestep degrres of freedom
	*/
	virtual VectorXd step(VectorXd& z, sim_state& state, VectorXd& f_ext, VectorXd& bc)
	{
		fast_cd_arap_static_precomp* sp = (fast_cd_arap_static_precomp*)this->sp;
		////needs to be updated every timestep
		assert(sp->AeqB.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		assert(f_ext.rows() == sp->BAB.rows() && "Force needs  to be be same dimensinoality as system we're solving");
		assert(z.rows() == sp->BAB.rows() && "intiial guess must be same dimensinoality as system we're solving");
		VectorXd z_next = z;

		VectorXd p = VectorXd::Zero(0);
		VectorXd x = Map<VectorXd>(params->X.data(), params->X.rows()* params->X.cols());
		
		((fast_cd_arap_dynamic_precomp*)dp)->precomp(z, p, state, f_ext, bc, *sp);
		z_next = ((fast_ik_arap_local_global_solver*)sol)->solve(z_next, *((fast_cd_sim_params*)params), *((fast_cd_arap_dynamic_precomp*)dp), *sp);
		return z_next;
	}

	/*
	Updates the equality constraints to new ones as given by matrix Aeq
	*/
	virtual void set_equality_constraint(SparseMatrix<double>& Aeq)
	{
		((fast_cd_sim_params*)params)->set_equality_constraint(Aeq);
		((fast_cd_arap_static_precomp*)sp)->set_equality_constraint(*((fast_cd_sim_params*)params));
		sol = new fast_ik_arap_local_global_solver(
			((fast_cd_arap_static_precomp*)sp)->BAB, ((fast_cd_arap_static_precomp*)sp)->AeqB,
			sol->p);
	}


};