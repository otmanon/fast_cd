#pragma once
#include "corot_precomp.h"
#include "fast_sim_params.h"
#include "repeat_for_each_entry.h"
#include "repeat_mat_for_each_entry.h"

#include "grouping_matrix_from_clusters.h"
#include <igl/diag.h>

using namespace Eigen;
struct fast_corot_static_precomp : public corot_static_precomp
{
	//modes
	MatrixXd BCmuB;
	MatrixXd BClamB;

	VectorXd BCmux;
	VectorXd BClamx;

	MatrixXd BAB;


	MatrixXd BMB;

	//clusters
	SparseMatrix<double> G; // basic grouping matrix
	SparseMatrix<double> G1;
	SparseMatrix<double> Gw;
	VectorXd GKx;

	//clusters + modes mix!
	MatrixXd GKB;
	MatrixXd GVmuKB;
	MatrixXd GVlamKB;

	VectorXd cluster_vols;
	VectorXd cluster_lambda;
	VectorXd cluster_mu;

	fast_corot_static_precomp() {};

	fast_corot_static_precomp(fast_sim_params& p) : corot_static_precomp(p)
	{
		
		BCmuB = p.B.transpose() * Cmu *p.B;
		BClamB = p.B.transpose() * Clam * p.B;

		BCmux = p.B.transpose() * Cmux;
		BClamx = p.B.transpose() * Clamx;
		BAB = p.B.transpose() * A * p.B;

		BMB = p.B.transpose() * M * p.B;



		grouping_matrix_from_clusters(p.labels, G);
		
		//form cluster mass weighing matrix
		cluster_vols = G * tet_vols;
		VectorXd cluster_volsi = 1/ cluster_vols.array();
		SparseMatrix<double> Vci, Vt, Gm_small;

		igl::diag(cluster_volsi, Vci);
		igl::diag(tet_vols, Vt);
		
		Gm_small = Vci * G * Vt;
		cluster_lambda = Gm_small * p.lambda;
		cluster_mu = Gm_small * p.mu;
		Gw = repeat_mat_for_each_entry(Gm_small, 3, 3);
		

		G1 = repeat_mat_for_each_entry(G, 3, 3);
	
		
		GKB = Gw * K * p.B;
		GKx = Gw * Kx;
		GVmuKB = G1 * VmuK * p.B;
		GVlamKB = G1 * VlamK * p.B;
	 //Eigen::SparseMatrix<double> G, G_tmp, S_cols, S_rows, G_1, G_m;
		// grouping_matrix_from_clusters(sol_p.labels, G);
	}
};

struct fast_corot_dynamic_precomp :corot_dynamic_precomp
{
	VectorXd BMBy;

	void precomp(const VectorXd& z,
		const sim_state& st, const VectorXd& f_ext,
		const VectorXd& bc, const fast_corot_static_precomp& sp)
	{
		assert(z.rows() == sp.BAB.rows() && "d.o.f.s need to match precomputation dimesnionlaity");
		assert(f_ext.rows() == z.rows() && "Force needs to have same number of rows as D.O.Fs");
		assert(bc.cols() == 1 && "Need to have at least one column (even if no rows)");

		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;

		BMBy = sp.BMB * z_hist;

	}
};