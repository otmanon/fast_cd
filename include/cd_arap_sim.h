#pragma once
#include "deformation_gradient_from_u_prefactorized_matrices.h"

#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std; using namespace igl; using namespace Eigen;

struct cd_arap_sim
{
public :

	cd_sim_params params;

	local_global_solver solver;

	cd_arap_static_precomp sp;

	cd_arap_dynamic_precomp dp;

	void init(cd_sim_params& params, local_global_solver& solver, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		this->params = params;
		this->solver = solver;
		this->sp = sp;
		this->dp = dp;
	}

	void step(VectorXd& z, VectorXd& p, cd_sim_state state, VectorXd& f_ext, VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params.Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");
		dp.precomp(z, p, state, f_ext, bc, sp);
		z = solver.solve(z, params, dp, sp);
	}
};

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
	VectorXd Mx;

	//C*J
	SparseMatrix<double> CJ;
	VectorXd Cx;

	void precomp(cd_sim_params& p)
	{
		VectorXd x = Map<VectorXd>(p.X.data(), p.X.rows() * p.X.cols());
		igl::massmatrix(p.X, p.T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		
		VectorXd Inull;
		deformation_gradient_from_u_prefactorized_matrices(p.X, p.T, Inull, K, V);
		
		SparseMatrix<double> Mu;
		igl::diag(p.mu, Mu);
		
		C = K.transpose() * Mu * V * K;
		CJ = C * p.J;
		Cx = C * x;

		SparseMatrix<double> C_test;
		igl::cotmatrix(p.X, p.T, C_test);
	
		C_test = igl::repdiag(C_test, 3);
		cout << "Cotan Laplacian Check" << (C - C_test).cwiseAbs().sum() << endl;
		
		SparseMatrix<double> I = SparseMatrix<double>(C.rows(), C.cols());
		A = (p.do_inertia * p.invh2 * M + C + p.Q + (!p.do_inertia) * 1e-9 * I);
		
		KJ = K * p.J;
		Kx = Mu * K * x;

		VK = V * Mu * K;
		VKJ = V * Mu * KJ;
		VKJ = V * Mu * Kx;

		MJ = M * p.J;
		Mx = M * x;
	}
};

struct cd_arap_dynamic_precomp
{
	VectorXd My;

	VectorXd Mur;
	VectorXd Cur;

	VectorXd VKur;

	//boundary conditions
	VectorXd bc;

	// Body forces
	VectorXd f_ext;

	void precomp(VectorXd& z, VectorXd& p, 
		cd_sim_state& st, VectorXd& f_ext, 
		VectorXd& bc, cd_arap_static_precomp& sp)
	{
		this->f_ext = f_ext;
		this->bc = bc;

		VectorXd z_hist = 2.0 * st.z_curr - st.z_prev;
		VectorXd p_hist = 2.0 * st.p_curr - st.p_prev;
		
		Cur = sp.CJ * p - sp.Cx;
		Mur = sp.MJ * p - sp.Mx;
		VectorXd Mur_hist = sp.MJ * p_hist - sp.Mx;
		My = Mur_hist + sp.M * z_hist - Mur;

	}

};


struct local_global_solver
{
	local_global_solver_params p;
	VectorXd z;
	min_quad_with_fixed_data<double> data;
	void init(SparseMatrix<double>& A, SparseMatrix<double>& Aeq, local_global_solver_params& p)
	{
		this->p = p;
		VectorXi bI;
		min_quad_with_fixed_precompute(A, bI, Aeq, true, data);
	}

	VectorXd solve(VectorXd& z, cd_sim_params& params, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		VectorXd r = local_step(z, dp, sp);
		z = global_step(z, params, dp, sp, r);
	};

	VectorXd local_step(VectorXd& z, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp)
	{
		VectorXd f = sp.VK * z + dp.VKur + sp.VKx;
		int nt = f.rows() / 9;
		MatrixXd F_stack = Map<MatrixXd>(f.data(), nt * 3, 3);
		MatrixXd R = MatrixXd::Zero(F_stack.rows(), F_stack.cols());

		Matrix3d F, rot;
		for (int i = 0; i < nt; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		VectorXd r = Map<const VectorXd>(R.data(), R.rows() * R.cols());
		return r;
	}

	VectorXd global_step(VectorXd& z, cd_sim_params& params, cd_arap_dynamic_precomp& dp, cd_arap_static_precomp& sp, VectorXd& r)
	{
		VectorXd inertia_grad = -dp.My;
		VectorXd arap_grad = dp.Cur + sp.Cx - sp.VK.transpose() * r;

		VectorXd g = params.invh2 * params.do_inertia * inertia_grad + arap_grad;
		VectorXd Y;

		VectorXd z_next;
		min_quad_with_fixed_solve(data, g, Y, dp.bc, z_next);
		//arap = dpre.Lur + spre.Lx - spre.MK' * r;
		return z_next;
	}

};
struct local_global_solver_params
{

	int max_iters;
	bool to_convergence;
	double threshold;
};


struct cd_sim_params
{
	//geometry
	MatrixXd X;
	
	//tet indices
	MatrixXi T;

	//timestep
	double h;
	double invh2;
	
	// first lame parameter per tet
	VectorXd mu;

	//second lame parameter per tet
	VectorXd lambda;

	//Equality constraints given as input from user
	SparseMatrix<double> Aeq;

	// linear rig jacobian! #V*dim x #bones*dim*dim+1
	SparseMatrix<double> J;

	//whether to activate inertia or not. if false, will add a bit of regularization to quadratic term
	bool do_inertia;

	//Extra user defined quadratic term (like for mass springs)
	SparseMatrix<double> Q;

	void default_init(MatrixXd& X, MatrixXi& T, 
		SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia )
	{
		this->X = X;
		this->T = T;
		this->mu = mu* VectorXd::Ones(X.rows() * 3);
		this->lambda = VectorXd::Zero(X.rows() * 3);
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;

		Q.resize(X.rows() * 3, X.rows() * 3);
		Q.setZero();

		Aeq = J.transpose().eval(); // no D-matrix here!
	}
	
};


struct cd_sim_state
{
	VectorXd z_curr, z_prev, p_curr, p_prev;

	void init(VectorXd& z_curr, VectorXd& z_prev, VectorXd& p_curr, VectorXd& p_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = p_curr;
		this->p_prev = p_prev;
	}

	void update(VectorXd& z_next, VectorXd& p_next)
	{
		this->z_prev = z_curr;
		this->z_curr = z_next;
		this->p_prev = p_curr;
		this->p_curr = p_next;
	}
};