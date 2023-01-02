#pragma once
#include "cd_sim_params.h"
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "cd_sim_state.h"

#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/volume.h>
#include <igl/diag.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;

struct cd_arap_static_precomp
{
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


	VectorXd tet_vols;

	cd_arap_static_precomp() {};


	cd_arap_static_precomp(cd_sim_params& p)
	{
		VectorXd x = Map<VectorXd>(p.X.data(), p.X.rows() * p.X.cols());
		igl::massmatrix(p.X, p.T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		Eigen::VectorXd m = M.diagonal();
		VectorXd Inull;
		deformation_gradient_from_u_prefactorized_matrices(p.X, p.T, Inull, K, V);
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
			CJ  = SparseMatrix<double>(C.rows(), 0);
			CJ.setZero();
		}
		Cx = C * x;

		/*SparseMatrix<double> C_test;
		igl::cotmatrix(sol_p.X, sol_p.T, C_test);
		C_test = igl::repdiag(C_test, 3);
		cout << "Cotan Laplacian Check " << (C + C_test).cwiseAbs().sum() << endl;*/

		SparseMatrix<double> I = SparseMatrix<double>(C.rows(), C.cols());
		I.setIdentity();
		A = (p.do_inertia * p.invh2 * M + C +  (!p.do_inertia) * 1e-9 * I);


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

	}


	virtual bool read_from_cache(string& precomp_cache_dir)
	{
		return false;
	};

	virtual bool write_to_cache(string& precomp_cache_dir)
	{
		return false;
	}

};


struct cd_arap_dynamic_precomp
{
	VectorXd My;

	VectorXd Mur;
	VectorXd Cur;

	VectorXd VKur;
	VectorXd Kur;
	//boundary conditions
	VectorXd bc;

	// Body forces
	VectorXd f_ext;

	void precomp(const VectorXd& z, const VectorXd& p,
		const cd_sim_state& st, const VectorXd& f_ext,
		const VectorXd& bc, const cd_arap_static_precomp& sp)
	{
		assert(f_ext.rows() == z.rows() && "Force needs to have same number of rows as D.O.Fs");
		assert(bc.cols() == 1 && "Need to have at least one column (even if no rows)");

		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;
		VectorXd p_hist = 2.0 * st.p_curr - st.p_prev;

		Cur = sp.CJ * p - sp.Cx;
		Mur = sp.MJ * p - sp.Mx;
		VKur = sp.VKJ * p - sp.VKx;
		Kur = sp.KJ * p - sp.Kx;

		VectorXd Mp_hist = (sp.MJ * p_hist);
		Eigen::VectorXd rig_momentum_terms = Mp_hist - sp.MJ * p;
		My = sp.M * z_hist + rig_momentum_terms;

	}

	//no rig params
	void precomp(const VectorXd& z,
		const cd_sim_state& st,const  VectorXd& f_ext,
		const VectorXd& bc, const cd_arap_static_precomp& sp)
	{
		assert(f_ext.rows() == z.rows() && "Force needs to have same number of rows as D.O.Fs");
		assert(bc.cols() == 1 && "Need to have at least one column (even if no rows)");

		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;

		Cur = VectorXd::Zero(z.rows());
		Mur = VectorXd::Zero(z.rows());
		VKur = VectorXd::Zero(sp.V.rows());
		Kur = VectorXd::Zero(sp.V.rows());

		My = sp.M * z_hist;// -sp.Mx;

	}

};

