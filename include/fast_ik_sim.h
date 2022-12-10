#pragma once
#include "fast_cd_arap_sim.h"
#include "lbs_jacobian.h"
#include "selection_matrix.h"

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

		VectorXi bIx = bI, bIy = bI, bIz = bI;
		bIy = bIx.array() + V.rows();
		bIz = bIy.array() + V.rows();
		VectorXi I = igl::cat(1, bIx, igl::cat(1, bIy, bIz));
		int n = V.rows() * 3;
		SparseMatrix<double> S;
	
		selection_matrix(bI, n, S);
		//SparseMatrix<double> J(0, B.rows());
		S = S.transpose().eval();
		params = new fast_cd_sim_params(V, T,
			Bd, l, S, 1, 0, 1e-2, false, "custom");

		fast_cd_arap_static_precomp* fcd_sp = new fast_cd_arap_static_precomp(*((fast_cd_sim_params*)params));

		sp = fcd_sp; // fast_cd_arap_static_precomp(sim_params);
		dp = new fast_cd_arap_dynamic_precomp();

		sol = new fast_cd_arap_local_global_solver(fcd_sp->BAB, fcd_sp->AeqB, solver_params);

	};
};