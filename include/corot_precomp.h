#pragma once
#include "cd_sim_params.h"
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "sim_state.h"
#include "vector_gradient_operator.h"
#include "repeat_for_each_entry.h"

#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/volume.h>
#include <igl/diag.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;

struct corot_static_precomp
{
	//cotan laplacian with youngs modulus
	SparseMatrix<double>  Cmu;

	//cotan laplacian with inccompressibility
	SparseMatrix<double>  Clam;

	VectorXd Cmux;
	VectorXd Clamx;

	//mass matrix
	SparseMatrix<double>  M;

	// Gradient matrix
	SparseMatrix<double> K;

	// tetrahedron volume matrix! This is NOT Vertex positions
	SparseMatrix<double> V;

	//system matrix
	SparseMatrix<double> A;

	//K*J, component of deformation gradient that comes from rig
	VectorXd Kx;

	// VKJ mass weighed  per tet KJ
	SparseMatrix<double> VmuK;
	VectorXd VmuKx;

	SparseMatrix<double> VlamK;
	VectorXd VlamKx;

	//M*x
	VectorXd Mx;

	VectorXd tet_vols;

	corot_static_precomp() {};


	corot_static_precomp(corot_sim_params& p)
	{
		VectorXd x = Map<VectorXd>(p.X.data(), p.X.rows() * p.X.cols());
		igl::massmatrix(p.X, p.T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		Eigen::VectorXd m = M.diagonal();
		vector_gradient_operator(p.X, p.T, K);

		igl::volume(p.X, p.T, tet_vols);

		//should abstract this. 
		VectorXd exp_vols = repeat_for_each_entry(tet_vols, 3, 3);
		igl::diag(exp_vols, V);

		SparseMatrix<double> Mu, Lam;
		VectorXd exp_mu, exp_lam;
		exp_mu = repeat_for_each_entry(p.mu, 3, 3);
		exp_lam = repeat_for_each_entry(p.lambda, 3, 3);
		igl::diag(exp_mu, Mu);
		igl::diag(exp_lam, Lam);

		SparseMatrix<double> I = SparseMatrix<double>(M.rows(), M.cols());
		I.setIdentity();
		Cmu =  K.transpose() * Mu * V * K + (!p.do_inertia) * 1e-9 * I;
		Clam = K.transpose() * Lam * V * K;
		Cmux = Cmu * x;
		Clamx = Clam * x;


		A = p.do_inertia * p.invh2 * M + Cmu;
		Kx = K * x;
		VmuK = V * Mu * K;
		VlamK = V * Lam * K;
		VmuKx = V * Mu * Kx;
		VlamKx = V * Lam * Kx;


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


struct corot_dynamic_precomp
{
	VectorXd My;

	
	//boundary conditions
	VectorXd bc;

	// Body forces
	VectorXd f_ext;

	void precomp(const VectorXd& z,
		const sim_state& st, const VectorXd& f_ext,
		const VectorXd& bc, const corot_static_precomp& sp)
	{
		assert(f_ext.rows() == z.rows() && "Force needs to have same number of rows as D.O.Fs");
		assert(bc.cols() == 1 && "Need to have at least one column (even if no rows)");

		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;
		
		My = sp.M * z_hist;

	}



};

