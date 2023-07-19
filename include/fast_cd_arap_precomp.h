#pragma once
#include "cd_arap_precomp.h"
#include "fast_cd_arap_sim_params.h"
#include "grouping_matrix_from_clusters.h"
#include "interweaving_matrix.h"
#include "read_fast_cd_sim_static_precompute.h"
#include <filesystem>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/sum.h>
struct fast_cd_arap_static_precomp : cd_arap_static_precomp
{
	/////////////////// FULL SPACE PRECOMP //////////////////
	//cotan laplacian
	SparseMatrix<double>  C;

	//mass matrix
	SparseMatrix<double>  M;

	// Gradient matrix
	SparseMatrix<double> K;

	// tetrahedron volume matrix! This is NOT Vertex positions
	SparseMatrix<double> V;

	//system matrix
	SparseMatrix<double> A;

	//K*J, component of deformation gradient that comes from rig
	SparseMatrix<double> KJ;
	VectorXd Kx;

	SparseMatrix<double> VK;
	// VKJ mass weighed  per tet KJ
	SparseMatrix<double> VKJ;
	VectorXd VKx;

	//M*J
	SparseMatrix<double>  MJ;

	//J'MJ
	MatrixXd  JMJ;


	//M*x
	VectorXd Mx;
	// vector here but really this is a scalar... only a vector for eigen consistence
	VectorXd xMx;
	//J'*M*x
	VectorXd JMx;

	//C*J
	SparseMatrix<double> CJ;
	VectorXd Cx;

	// tet volumes
	VectorXd tet_vols;


	////////// SUBSPACE PRECOMP ///////////////////////
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


	fast_cd_arap_static_precomp(fast_cd_arap_sim_params& p) 
	{
		init(p);
	}	


	
	void init(fast_cd_arap_sim_params& p)
	{

		/////////////////// FULL SPACE PRECOMP ////////////////
		VectorXd x = Map<VectorXd>(p.X.data(), p.X.rows() * p.X.cols());
		igl::massmatrix(p.X, p.T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		Eigen::VectorXd m = M.diagonal();
		VectorXd Inull;
		deformation_jacobian(p.X, p.T, Inull, K, V);
		Eigen::VectorXd v = V.diagonal();

		igl::volume(p.X, p.T, tet_vols);

		SparseMatrix<double> Mu;
		igl::diag(p.mu, Mu);
		Mu = igl::repdiag(Mu, 9);

		C = K.transpose() * Mu * V * K;
		if (p.J.rows() > 0)
			CJ = C * p.J;
		else
		{
			CJ = SparseMatrix<double>(C.rows(), 0);
			CJ.setZero();
		}
		Cx = C * x;
		SparseMatrix<double> I = SparseMatrix<double>(C.rows(), C.cols());
		I.setIdentity();
		A = (p.do_inertia * p.invh2 * M + C + (!p.do_inertia) * 1e-9 * I);


		if (p.J.rows() > 0)
			KJ = K * p.J;
		else
		{
			KJ = SparseMatrix<double>(K.rows(), 0);
			KJ.setZero();
		}
		Kx = K * x;

		VK = V * Mu * K;
		VKJ = V * Mu * KJ;

		VKx = V * Mu * Kx;

		if (p.J.rows() > 0)
			MJ = M * p.J;
		else
		{
			MJ = SparseMatrix<double>(M.rows(), 0);
			MJ.setZero();
		}
		Mx = M * x;

		JMJ = p.J.transpose() * M * p.J;
		JMx = p.J.transpose() * M * x;
		xMx = x.transpose() * M * x;
		
		///////////////////////// SUBSPACE PRECOMP ////////////////////////
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

	void set_equality_constraint(fast_cd_arap_sim_params& p)
	{
		AeqB = p.Aeq * p.B;
	}


	bool read_from_cache(string& precomp_cache_dir)
	{
		bool t = true;
		t = igl::readDMAT(precomp_cache_dir + "/BCB.DMAT", BCB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BMB.DMAT", BMB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BAB.DMAT", BAB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/AeqB.DMAT", AeqB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/GmKB.DMAT", GmKB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/GmKJ.DMAT", GmKJ) && t;
		t = igl::readDMAT(precomp_cache_dir + "/GmKx.DMAT", GmKx) && t;
		t = igl::readDMAT(precomp_cache_dir + "/G1VKB.DMAT", G1VKB) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BMJ.DMAT", BMJ) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BMx.DMAT", BMx) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BCJ.DMAT", BCJ) && t;
		t = igl::readDMAT(precomp_cache_dir + "/BCx.DMAT", BCx) && t;
		t = igl::readDMAT(precomp_cache_dir + "/JMx.DMAT", JMx) && t;
		t = igl::readDMAT(precomp_cache_dir + "/xMx.DMAT", xMx) && t;
		t = igl::readDMAT(precomp_cache_dir + "/JMJ.DMAT", JMJ) && t;
	
		if (!t)
			printf(" precomp cache dir %s, is either corrupt or outdated. please construct fast_cd_arap_sim differently \n");
		return t;
	}

	bool write_to_cache(string& precomp_cache_dir)
	{
		return false;
		/*namespace fs = std::filesystem;
		if (!fs::exists(fs::path(precomp_cache_dir)))
		{
			fs::create_directories(fs::path(precomp_cache_dir));
		}
	
		printf("Saving fast_cd_arap matrix precomputations to cache dir %s ... \n", precomp_cache_dir.c_str());
		bool t = true;
		t = igl::writeDMAT(precomp_cache_dir + "/BCB.DMAT", BCB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BMB.DMAT", BMB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BAB.DMAT", BAB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/AeqB.DMAT", AeqB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/GmKB.DMAT", GmKB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/GmKJ.DMAT", GmKJ) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/GmKx.DMAT", GmKx) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/G1VKB.DMAT", G1VKB) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BMJ.DMAT", BMJ) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BMx.DMAT", BMx) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BCJ.DMAT", BCJ) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/BCx.DMAT", BCx) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/JMx.DMAT", JMx) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/xMx.DMAT", xMx) && t;
		t = igl::writeDMAT(precomp_cache_dir + "/JMJ.DMAT", JMJ) && t;
		if (!t)
		{
			printf("Could not save fast_cd_arap precomp cache!");
		}
		return t;*/
	}
};

struct fast_cd_arap_dynamic_precomp : cd_arap_dynamic_precomp
{

	VectorXd BMy;

	VectorXd BMur;
	VectorXd BCur;

	VectorXd GmKur;

	void precomp(const VectorXd& z, const VectorXd& p,
		const cd_sim_state& st, const VectorXd& f_ext,
		const VectorXd& bc, const fast_cd_arap_static_precomp& sp)
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