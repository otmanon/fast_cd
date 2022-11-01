#pragma once
#include "cd_arap_precomp.h"
#include "fast_cd_sim_params.h"
#include "grouping_matrix_from_clusters.h"
#include "interweaving_matrix.h"
#include "read_fast_cd_sim_static_precompute.h"
#include <filesystem>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/sum.h>
struct fast_cd_arap_static_precomp : cd_arap_static_precomp
{
	//reduced cotan laplacian (with heterogeneity inside)
	MatrixXd BCB;

	//reduced mass matrix (sometimes will be identity if B'MB= I is enforce in eigs)
	MatrixXd BMB;

	//reduced system matrix
	MatrixXd BAB;

	//reduced equality constraint
	MatrixXd AeqB;

	

	//mass, stiffness weighed reduced  gradient operator
	MatrixXd GmKB;
	MatrixXd GmKJ;
	VectorXd GmKx;
	
	MatrixXd G1VKB;


	MatrixXd BMJ;
	VectorXd BMx;

	MatrixXd BCJ;
	VectorXd BCx;
	
	fast_cd_arap_static_precomp() {};


	fast_cd_arap_static_precomp(fast_cd_sim_params& p) : cd_arap_static_precomp(p)
	{
		init(p);
	}	


	
	void init(fast_cd_sim_params& p)
	{
		BCB = p.B.transpose() * C * p.B;
		BMB = p.B.transpose() * M * p.B;
		BAB = p.B.transpose() * A * p.B;
		AeqB = p.Aeq * p.B;

		//going from per cluster rotations to forces
		GmKB = K * p.B;
		GmKJ = KJ;
		GmKx = Kx;

		//going from displacements to per-cluster deformation gradients.
		G1VKB = VK * p.B;

		BMJ = p.B.transpose() * MJ;
		BMx = p.B.transpose() * Mx;

		BCJ = p.B.transpose() * CJ;
		BCx = p.B.transpose() * Cx;

		// if we have tet clusters, precompute corresponding matrices
		if (p.labels.maxCoeff() < p.T.rows() - 1 && p.labels.rows() > 0)
		{
			Eigen::SparseMatrix<double> G, G_tmp, S_cols, S_rows, G_1, G_m;
			grouping_matrix_from_clusters(p.labels, G);


			G_tmp = igl::repdiag(G, 3);
			interweaving_matrix(G.cols(), 3, S_cols);
			interweaving_matrix(G.rows(), 3, S_rows);
			G_tmp = S_rows.transpose() * G_tmp * S_cols;
			G_tmp = igl::repdiag(G_tmp, 3);
			G_1 = G_tmp;
			Eigen::VectorXd cluster_mass;

			G_tmp = G_1 * V;
			igl::sum(G_tmp, 2, cluster_mass);

			G_m = G_1 * V;
			igl::sum(G_m, 2, cluster_mass);


			cluster_mass.array() = 1.0 / cluster_mass.array();
			G_m = cluster_mass.asDiagonal() * G_m;

			//one is mass weighed, used to average per tet quantities to per cluster ones
			GmKB = G_m * GmKB;
			GmKJ = G_m * KJ;
			GmKx = G_m * Kx;
			//this one is an indicator, used to map per cluster quantities to each
			G1VKB = G_1 * G1VKB;
		}
	}
};

struct fast_cd_arap_dynamic_precomp : cd_arap_dynamic_precomp
{

	VectorXd BMy;

	VectorXd BMur;
	VectorXd BCur;

	VectorXd GmKur;

	void precomp(VectorXd& z, VectorXd& p,
		cd_sim_state& st, VectorXd& f_ext,
		VectorXd& bc, fast_cd_arap_static_precomp& sp)
	{
		assert(f_ext.rows() == z.rows() && "Force needs to have same number of rows as D.O.Fs");
		assert(bc.cols() == 1 && "Need to have at least one column (even if no rows)");
		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;
		VectorXd p_hist = 2.0 * st.p_curr - st.p_prev;

		BCur = sp.BCJ * p - sp.BCx;
		BMur = sp.BMJ * p - sp.BMx;
		GmKur = sp.GmKJ * p - sp.GmKx;

		VectorXd BMp_hist = (sp.BMJ * p_hist);
		Eigen::VectorXd rig_momentum_terms = BMp_hist - sp.BMJ * p;
		BMy = sp.BMB * z_hist + rig_momentum_terms;
	}
};