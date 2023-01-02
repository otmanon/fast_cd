#pragma once
#include "fast_corot_precomp.h"
#include "fast_cd_corot_sim_params.h"

struct fast_cd_corot_static_precomp : public fast_corot_static_precomp
{
	MatrixXd BMJ;
	MatrixXd BCmuJ;
	SparseMatrix<double> GKJ;
	
	fast_cd_corot_static_precomp(fast_cd_sim_params& p) :  fast_corot_static_precomp(p)
	{
		BMJ = p.B.transpose() * M* p.J;
		BCmuJ = p.B.transpose() * Cmu * p.J;
		GKJ = Gw * K * p.J;
	}

};

struct fast_cd_corot_dynamic_precomp : public fast_corot_dynamic_precomp
{
	VectorXd BMy;
	void precomp(const VectorXd& z, const VectorXd& p,
		const cd_sim_state& st, const VectorXd& f_ext,
		const VectorXd& bc, const fast_cd_corot_static_precomp& sp)
	{

		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;
		VectorXd p_hist = 2.0 * st.p_curr - st.p_prev;

		
		VectorXd BMp_hist = (sp.BMJ * p_hist);
		Eigen::VectorXd rig_momentum_terms = BMp_hist - sp.BMJ * p;
		BMy = sp.BMB * z_hist + rig_momentum_terms;
	}

};