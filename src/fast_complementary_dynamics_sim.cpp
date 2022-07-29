#include "fast_complementary_dynamics_sim.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/average_onto_faces.h>
#include <igl/polar_svd3x3.h>
#include <igl/slice.h>
#include <igl/cat.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <iostream>



#include "arap_hessian.h"
#include "trace_matrix_operator.h"
#include "covariance_scatter_matrix.h"
#include "deformation_gradient.h"
#include "interweaving_matrix.h"
#include "grouping_matrix_from_clusters.h"
#include "fast_complementary_dynamics_constrained_gevp.h"
#include "complementary_equality_constraint.h"
#include "compute_modes_matlab.h"
#include "compute_modes_spectra.h"
#include "compute_clusters_igl.h"
#include "compute_clusters_matlab.h"
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include <igl/sum.h>
#include "igl/count.h"
#include <igl/colon.h>
#include <igl/volume.h>
#include <igl/get_seconds.h>
#include <split_components.h>
#include <filesystem>

#include <igl/invert_diag.h>
#include "DenseQuasiNewtonSolver.h"
#include "QuasiNewtonSolver.h"

#include "per_vertex_rotation.h"

FastCDSim::FastCDSim( Eigen::MatrixXd& X,  Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& J, Eigen::MatrixXd& W, double ym, double pr,
	double dt, int num_modes, int num_clusters, std::string modes_file_dir, std::string clusters_file_dir, 
	bool do_reduction, bool do_clustering, int num_modal_features, bool do_inertia,
	double momentum_leaking_dt, std::string metric_type, double metric_alpha, double beta, double gamma)
{
	this->X = X;
	this->T = T;
	this->J = J;
	this->num_modes = num_modes;
	this->num_clusters = num_clusters;
	this->do_clustering = do_clustering;
	this->do_reduction = do_reduction;
	this->modes_file_dir = modes_file_dir;
	this->clusters_file_dir = clusters_file_dir;
	this->dt = dt;
	this->num_modal_features = num_modal_features;
	this->do_inertia = do_inertia;
	this->momentum_leaking_dt = momentum_leaking_dt;
	this->metric_type = metric_type;
	this->metric_alpha = metric_alpha;
	this->beta = beta;
	this->gamma = gamma;
	this->W = W;
	stiffness = ym / (2.0 * (1.0 + pr));
	incompressibility = ym* pr / ((1.0 + pr) * (1.0 - 2.0 * pr));

	constraints.use_constraints = false;

	Eigen::VectorXi _n; //throwaway var
	igl::boundary_facets(T, F, FiT, _n);
	//get list of exterior vertex indices
	igl::unique(F, ext_ind);
	
	reduced_newton_solver = new DenseQuasiNewtonSolver();
	full_newton_solver = new QuasiNewtonSolver();

	constraints.reduced_newton_solver = new DenseQuasiNewtonSolver();
	constraints.full_newton_solver = new QuasiNewtonSolver();
	update_compelementary_constraint(J, modes_file_dir, clusters_file_dir);
}

Eigen::VectorXd FastCDSim::reduced_step_with_equality_constriants(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& bc, bool to_convergence, int max_iters, double convergence_threshold)
{

	if (!constraints.use_constraints)
	{
		std::cout << "simulator is configured to NOT use equality constraints... please call reduced_step() instead" << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	if (!do_reduction)
	{
		std::cout << "simulator not configured for calling reduced_step(). Please switch simulator to reduced mode with sim.switch_reduced(true)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}

	Eigen::VectorXd r, ur, BMy;
	//have a valid rig
	BMy = rmp.BMB * (2.0 * z_curr - z_prev) +  gamma*rmp.BMJ * (2.0 * p_curr - p_prev);            //u_curr and u_prev are total displacements, eg u_prev = uc_prev + ur_prev;
//	dm.r = J * p_next;
	dm.BMy_tilde =  gamma* rmp.BMJ * p_next - BMy ;

//	printf("r: %r\n", dm.r.norm());//
	dm.GMKur = fmp.GMKJ * p_next - fmp.GMKX;
	dm.BKMKur = rmp.BKMKJ * p_next - rmp.BKMKX;

//	dm.BKMTIKur = rmp.BKMTIKJ * p_next - rmp.BKMTIKX;
	dm.GmKur = fmp.GmKJ * p_next - fmp.GmKX;

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> f;
	reduced_qnewton_energy_grad(f, grad_f);
	Eigen::VectorXd z_next;
	bool do_line_search = incompressibility < 1e-8 ? false : true;
	z_next = constraints.reduced_newton_solver->solve_with_equality_constraints(z_curr, f, grad_f, bc, do_line_search, to_convergence, max_iters, convergence_threshold);
	return z_next;

}//


Eigen::VectorXd FastCDSim::reduced_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, bool to_convergence, int max_iters, double convergence_threshold)
{
	if (constraints.use_constraints)
	{
		std::cout << "simulator is configured to use equality constraints... please call reduced_step_with_equality_constraints() instead" << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	if (!do_reduction)
	{
		std::cout << "simulator not configured for calling reduced_step(). Please switch simulator to reduced mode with sim.switch_reduced(true)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	Eigen::VectorXd r, ur, BMy;
	//have a valid rig
	BMy = rmp.BMB *  (2.0 * z_curr - z_prev) + gamma*rmp.BMJ * (2.0*p_curr - p_prev);            //u_curr and u_prev are total displacements, eg u_prev = uc_prev + ur_prev;
//	dm.r = J * p_next;g
	dm.BMy_tilde = gamma* rmp.BMJ * p_next - BMy;

//	BMy = rmp.BTMB * (2.0 * z_curr - z_prev);
//	dm.BMy_tilde =  - BMy;
	dm.GMKur = fmp.GMKJ * p_next  - fmp.GMKX;
	dm.BKMKur = rmp.BKMKJ * p_next  - rmp.BKMKX;
	
//	dm.BKMTIKur = rmp.BKMTIKJ * p_next - rmp.BKMTIKX;
	dm.GmKur = fmp.GmKJ * p_next  - fmp.GmKX;

	Eigen::VectorXd uc_bc;
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> f;
	reduced_qnewton_energy_grad(f, grad_f);

	Eigen::VectorXd z_next;
	bool do_line_search = incompressibility < 1e-8 ? false : true;
	z_next = reduced_newton_solver->solve(z_curr, f, grad_f, do_line_search, to_convergence, max_iters, convergence_threshold);

	return z_next;
}

Eigen::VectorXd FastCDSim::full_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev, bool to_convergence, int max_iters, double convergence_threshold)
{
	if (do_reduction)
	{
		std::cout << "simulator not configured for calling full_step. Please switch simulator to full mode with sim.switch_reduced(false)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	Eigen::VectorXd My, r;
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	r = J * p_next;
	//compute dynamic matrix at the start of the timestep
	dm.r = r;//TODO: delete this once we're done playing around with it... unnecessary
	//precompute_with_constraints AJ and MJ for more speed
	My = fmp.M * (2.0 * uc_curr - uc_prev) + gamma*fmp.MJ * (2.0 * p_curr - p_prev);

	dm.My_tilde = gamma * fmp.M * dm.r -  My  ;// +2 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	dm.GmKur = fmp.GmKJ * p_next - fmp.GmKX;
//	dm.y_tilde = y - dm.r;// +2.0 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	full_qnewton_energy_grad(energy, grad_f);

	Eigen::VectorXd uc_bc;
	Eigen::VectorXd uc_next;
	
	Eigen::VectorXd zero = Eigen::VectorXd::Zero(J.cols());

	bool do_line_search =  incompressibility < 1e-8 ? false : true;
	uc_next = full_newton_solver->solve_with_equality_constraints(uc_curr, energy,  grad_f, zero, do_line_search, to_convergence, max_iters, convergence_threshold);
	
	return uc_next;
}

Eigen::VectorXd FastCDSim::full_step_with_equality_constraints(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev,
	const Eigen::VectorXd& bc, bool to_convergence, int max_iters, double convergence_threshold)
{
	if (do_reduction)
	{
		std::cout << "simulator not configured for calling full_step. Please switch simulator to full mode with sim.switch_reduced(false)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	Eigen::VectorXd My, r;
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	r = J * p_next;
	//compute dynamic matrix at the start of the timestep
	dm.r = r;//TODO: delete this once we're done playing around with it... unnecessary
	//precompute_with_constraints AJ and MJ for more speed
	My = fmp.M * (2.0 * uc_curr - uc_prev) + gamma*fmp.MJ * (2.0 * p_curr - p_prev);
	dm.My_tilde = gamma *fmp.M * dm.r - My;// +2 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	
	dm.GmKur = fmp.GmKJ * p_next - fmp.GmKX;
	//	dm.y_tilde = y - dm.r;// +2.0 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	full_qnewton_energy_grad(energy, grad_f);

	Eigen::VectorXd uc_bc;
	Eigen::VectorXd uc_next;
	Eigen::VectorXd zero = Eigen::VectorXd::Zero(J.cols());
	Eigen::VectorXd b = igl::cat(1, zero, bc);
	bool do_line_search = incompressibility < 1e-8 ? false : true;
	uc_next = constraints.full_newton_solver->solve_with_equality_constraints(uc_curr, energy, grad_f, b, do_line_search, to_convergence, max_iters, convergence_threshold);
	return uc_next;
}

void FastCDSim::update_compelementary_constraint(const Eigen::SparseMatrix<double>& J, std::string new_modes_dir, std::string new_clusters_dir)
{
	this->J = J;
	if (this->J.cols() == 0)
	{
		this->J.resize(X.rows() * 3, 0);
		this->J.setZero();
	}
	update_modes_cache_dir(new_modes_dir);
	update_clusters_cache_dir(new_clusters_dir);

	init_modes(num_modes);
	if (do_clustering)
		init_clusters(num_clusters, num_modal_features);
	else
		igl::colon(0, T.rows() - 1, labels);
	num_clusters = labels.maxCoeff()+1;
	init_system_matrices();
	precompute_solvers();
}

void FastCDSim::precompute_solvers()
{
	if (!constraints.use_constraints)
	{
		if (do_reduction)
			reduced_newton_solver->precompute(rmp.BAB);
		else
		{
			full_newton_solver->precompute_with_equality_constraints(fmp.A, fmp.Aeq);
		}	

	}
	if(constraints.use_constraints)
	{
		if (do_reduction)
			constraints.reduced_newton_solver->precompute_with_equality_constraints(rmp.BAB, constraints.SB);
		else
		{
			Eigen::SparseMatrix<double> Heq;
			igl::cat(1, fmp.Aeq, constraints.S, Heq);

			constraints.full_newton_solver->precompute_with_equality_constraints(fmp.A, Heq);
		}
	}
}


void FastCDSim::reduced_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
	f = [&](const Eigen::VectorXd z)
	{
		
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat = (fmp.GmH + rmp.GmKB * z + beta*dm.GmKur);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::VectorXd volume_energy_tet = Eigen::VectorXd::Zero(l);
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F;
		//fit rotations
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;

			//calculate (tr R^T F -I)R
			double tr = (rot.transpose() * F).diagonal().sum() - 3;
			//tr = tr * tr;
			volume_energy_tet(i) = 0.5 * fmp.Vol_c.coeff(i, i) * tr * tr;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

		double inertia = 0.5 / (dt * dt) * z.transpose() * rmp.BMB * z;
		inertia += 1.0 / (dt * dt) * z.transpose() * dm.BMy_tilde;
			//- 2.0 * z.transpose() * dm.BMr;
			//_y.transpose() * fmp.M * u_y;

//		Eigen::VectorXd f_r = fmp.H + fmp.K * z + dm.Kur - fmp.G_1.transpose() * R_flat;

	//	double bending = 0.5 * f_r.transpose() * fmp.Vol_exp * f_r;

		double bending =0.5 * z.transpose() * ( rmp.BKMH + rmp.BKMKB * z + beta*dm.BKMKur - rmp.BKMG * R_flat);
		//Eigen::VectorXd tr = fmp.traceMat * f_r;
		//double volume = 0.5 * tr.transpose() * fmp.Vol * tr;//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		double elastic = stiffness * bending + incompressibility * (volume_energy_tet.sum());

		double e = elastic + do_inertia * inertia;
		return e;

	};
	grad_f = [&](const Eigen::VectorXd z)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat = (fmp.GmH + rmp.GmKB * z + beta*dm.GmKur);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F;


		//fit rotations
		for (int i = 0; i < l; i++)
		{
			F = (F_stack.block(3 * i, 0, 3, 3));
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;// .transpose();//

	
			//calculate (tr R^T F -I)R
			const double tr = (rot.transpose() * F).diagonal().sum() - 3;
			trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
		}
	
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
		//	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;
		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB * z + dm.BMy_tilde);
		//Eigen::VectorXd bending_grad =  rmp.BKMG * (FV_flat - R_flat);
		
		Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB * z - rmp.BKMG * R_flat + beta*dm.BKMKur;
		//Eigen::VectorXd volume_grad = rmp.BKMTIH +  rmp.BKMTIKB* z + dm.BKMTIKur  - rmp.BKMTIG * R_flat; //sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		Eigen::VectorXd volume_grad = rmp.BKMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad + incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + do_inertia * inertia_grad;
	    
	//	printf("inertia: %g, arap: % g  volume: %g", inertia_grad.norm(), bending_grad.norm(), volume_grad.norm());
		return g;
	};
}

void FastCDSim::full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{

	f = [&](const Eigen::VectorXd uc)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd u = uc +  dm.r - Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

		double bending, volume, inertia;
		energy(u, dm.u_curr, dm.u_prev, bending, volume, inertia);
		double e = bending + volume + do_inertia * inertia;
		return e;
	};

	grad_f = [&](const Eigen::VectorXd uc)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd ur = dm.r - Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
		Eigen::VectorXd u = uc + ur;
		Eigen::VectorXd FV_flat = (fmp.GmH + fmp.GmK * (uc + beta*ur));//fmp.GMH + fmp.GMK*u;//sm.G_exp * sm.FM * (sm.H + sm.K * u);
        Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F; 
		Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		//fit rotations
		for (int i = 0; i < l; i++)
		{
			F = (F_stack.block(3 * i, 0, 3, 3));
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;

			const double tr = (rot.transpose() * F).diagonal().sum() - 3;
			trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());

		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (fmp.M *uc + dm.My_tilde);
		//Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB * z - dm.BMy_tilde);
		Eigen::VectorXd bending_grad = fmp.KMH + fmp.KMK * (uc + beta*ur) - fmp.KMG * R_flat;
	//	Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB * z + dm.BKMKur - rmp.BKMG * R_flat;

		Eigen::VectorXd volume_grad = fmp.KMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad + incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + do_inertia * inertia_grad;


	//const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
	//const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
	////	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;
	//Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB * z + dm.BMy_tilde);
	////Eigen::VectorXd bending_grad =  rmp.BKMG * (FV_flat - R_flat);
	//
	//Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB * z - rmp.BKMG * R_flat + beta * dm.BKMKur;
	////Eigen::VectorXd volume_grad = rmp.BKMTIH +  rmp.BKMTIKB* z + dm.BKMTIKur  - rmp.BKMTIG * R_flat; //sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//Eigen::VectorXd volume_grad = rmp.BKMG * gb_flat;
	//Eigen::VectorXd elastic_grad = stiffness * bending_grad + incompressibility * volume_grad;
	//Eigen::VectorXd g = elastic_grad + do_inertia * inertia_grad;
		return g;
	};
}

/*
Makes sure the simulation enforces positional constraints on mesh indices bi. bc is a list of desired vertex positions which will pe 
enforce as hard constraints in the optimizaiton
*/
void FastCDSim::make_constraints(Eigen::SparseMatrix<double> S)
{
	constraints.S = S;
	constraints.SB = constraints.S * B;

	
	//FullPivLU<MatrixXd> lu_decomp3(B.transpose() * B);
	//std::cout << "rank of B : " << lu_decomp3.rank() << std::endl;
	//std::cout << "number of row of B : " << B.transpose().cols() << std::endl;
	//
	//
	//Eigen::MatrixXd ST = S.transpose();
	//FullPivLU<MatrixXd> lu_decomp(S.transpose() * S);
	//std::cout << "rank of J : " << lu_decomp.rank() << std::endl;
	//std::cout << "number of row of J : " << S.transpose().rows() << std::endl;
	//
	//
	//FullPivLU<MatrixXd> lu_decomp2(constraints.SB.transpose() * constraints.SB);
	//std::cout << "rank of JB : " << lu_decomp2.rank() << std::endl;
	//std::cout << "number of row of JB : " << constraints.SB.rows() << std::endl;
	constraints.SJ = constraints.S * J;
	constraints.SX = constraints.S * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows()*X.cols());
	constraints.use_constraints = true;
	precompute_solvers();
}


void FastCDSim::init_modes(int num_modes)
{
	namespace fs = std::filesystem;

	std::string B_file_path = modes_file_dir + "B.DMAT";
	std::string L_file_path = modes_file_dir + "L.DMAT";
	bool found_modes = igl::readDMAT(B_file_path, B_full);
	bool found_evals = igl::readDMAT(L_file_path, L_full);
	bool enough_modes = false;
	
	if (found_modes && found_evals)
	{
		printf("Found %i / %i cached modes at %s", B_full.cols(), num_modes, B_file_path.c_str());
		enough_modes = num_modes <= B_full.cols();

	}
	if (!enough_modes)
	{
		printf("Could not find cached modes at %s, computing %i modes from scratch...\n", B_file_path.c_str(), num_modes);
		//build constrained eigenvalue problem
		Eigen::SparseMatrix<double> H;
		Eigen::SparseMatrix<double> M;
		fast_complementary_dynamics_constrained_geivp(X, T, J, H, M, momentum_leaking_dt, metric_type, metric_alpha);

		//compute modes and load into S_full and B_full
		Eigen::MatrixXd modes;
		#ifdef FAST_CD_USE_MATLAB
		//	compute_modes_spectra(H, M, num_modes, modes, L_full);
				compute_modes_matlab(H, M, num_modes, modes, L_full);
		#else
				compute_modes_spectra(H, M, num_modes, modes, L_full);
		#endif
		B_full = modes.block(0, 0, 3*X.rows(), num_modes);
		if (!fs::exists(fs::path(B_file_path).parent_path()))
		{
			fs::create_directories(fs::path(B_file_path).parent_path());
		}
		printf("Saving new modes at %s ... \n", B_file_path.c_str());
		//writeDMAt
		igl::writeDMAT(B_file_path, B_full, false);
		Eigen::MatrixXd L_full_mat = Eigen::Map<Eigen::MatrixXd>(L_full.data(), num_modes, 1); //annoying
		igl::writeDMAT(L_file_path, L_full_mat, false);
	}
	L = L_full.topRows(num_modes);
	B = B_full.block(0, 0, 3*X.rows(), num_modes);
	fmp.B_0 = B;
}

void FastCDSim::init_clusters(int num_clusters, int num_feature_modes)
{

	namespace fs = std::filesystem;

	const int l = do_clustering ? num_clusters : T.rows();
	std::string labels_file_path = clusters_file_dir + "labels_" + std::to_string(l) + "_features_" + std::to_string(num_feature_modes) + ".DMAT";
	std::string centroids_file_path = clusters_file_dir + "centroids_" + std::to_string(l) + "_features_" + std::to_string(num_feature_modes) + ".DMAT";

	bool found_clusters = igl::readDMAT(labels_file_path, labels);
	if (!found_clusters)
	{
		if (!(l == T.rows()))
		{
			printf("Could not find cached clusters at %s,  clustering...\n", labels_file_path.c_str());
			double t_start = igl::get_seconds();

			Eigen::MatrixXd C;
		#ifdef FAST_CD_USE_MATLAB
			//compute_clusters_matlab(T, B_full, L_full, num_clusters, num_feature_modes, labels, C);
			compute_clusters_igl(T, B_full, L_full, num_clusters, num_feature_modes, labels, C);
		#else
			compute_clusters_igl(T, B_full, L_full, num_clusters, num_feature_modes, labels, C);
		#endif
			printf("Done clustering! took %g seconds \n", igl::get_seconds() - t_start);

			if (!fs::exists(fs::path(labels_file_path).parent_path()))
			{
				fs::create_directories(fs::path(labels_file_path).parent_path());
			}
			printf("Saving at %s...\n", labels_file_path.c_str());
			Eigen::MatrixXi labels_mat = Eigen::Map<Eigen::MatrixXi>(labels.data(), labels.rows(), 1);
			igl::writeDMAT(labels_file_path, labels_mat, false);
			igl::writeDMAT(centroids_file_path, C, false);
			labels = labels_mat.col(0);
			
		}
		else {
			igl::colon(0, T.rows() - 1, labels);
		}
	}
	else
	{
		std::cout << "Found " << std::to_string(num_clusters) << " Cached Clusters at " << labels_file_path << std::endl;
		//build up grouping matrix dimesnion by kron Id
	}
	Eigen::SparseMatrix<double> G_small, G_large, S_cols, S_rows;
	grouping_matrix_from_clusters(labels, fmp.G);

}

void FastCDSim::init_full_system_matrices()
{
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	//cotan laplacian
	igl::cotmatrix(X, T, fmp.C);

	fmp.C = -igl::repdiag(fmp.C, 3);

	//mass matrix
	igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, fmp.M);
	fmp.M = igl::repdiag(fmp.M, 3);

	//constructs A matrix, and does preomputation if we will use it directly
	if (do_inertia)
		fmp.A = stiffness * fmp.C + (1.0 / dt) * (1.0 / dt) * fmp.M;
	else
	{
		Eigen::SparseMatrix<double> I(fmp.M.rows(), fmp.M.rows());
		I.setIdentity();
		fmp.A = stiffness * fmp.C + 1e-10 * I; // tiny regularization
	}


	complementary_equality_constraint(X, T, J, fmp.Aeq, momentum_leaking_dt, metric_type, metric_alpha);

	fmp.MJ = fmp.M * J;
	fmp.MX = fmp.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	//flattened_deformation_gradient matrices. 
	deformation_gradient_from_u_prefactorized_matrices(X, T, fmp.H, fmp.K, fmp.Vol_exp);

	Eigen::VectorXd vol;
	igl::volume(X, T, vol);
	fmp.Vol = vol.asDiagonal();
	trace_matrix_operator(T.rows(), 3, fmp.traceMat);
	
	
	//Need these for bending energy
	Eigen::SparseMatrix<double> KT = fmp.K.transpose();
	Eigen::SparseMatrix<double> KM = KT * fmp.Vol_exp;
	fmp.KMH = KM * fmp.H;
	fmp.KMK = KM * fmp.K;				//this should be exactly equal to cotan laplacian//
	
	Eigen::SparseMatrix<double> KMTT = KM * fmp.traceMat.transpose() * fmp.traceMat;
	fmp.KMTIH  = KMTT *  fmp.H;
	fmp.KMTIK =  KMTT * fmp.K;

	//should do clustering matrices here.

	fmp.KMKJ = fmp.KMK * J; 
	fmp.KMKX = fmp.KMK * x;


	fmp.JMJ = J.transpose() * fmp.M * J;
		


	//build H
	// build K
}

void FastCDSim::init_reduction_matrices()
 {
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols()); //a few of these need rest space positions
	//general system matrices
	Eigen::MatrixXd BT = B.transpose();
	Eigen::SparseMatrix<double> KT = fmp.K.transpose();
	rmp.BA = BT * fmp.A;
	rmp.BAB = rmp.BA * B;
	rmp.BAJ = rmp.BA * J;

	rmp.BM = BT * fmp.M;
	rmp.BMB = rmp.BM * B;
	rmp.BMJ = rmp.BM * J;

	//rmp.BMX = rmp.BM * x;

	//need these for bending energy
	rmp.BKMK = BT * fmp.KMK;
	rmp.BKMKB = rmp.BKMK * B;
	rmp.BKMH = BT * fmp.KMH;

	//rmp.GMKB = fmp.GMK * B;
	rmp.GmKB = fmp.GmK * B;
	rmp.BKMG = BT* fmp.KMG;

//	rmp.BKMGm = BT* (KT) * (fmp.Vol_exp * fmp.G_m.transpose());
	rmp.BKMKJ = rmp.BKMK * J;
	rmp.BKMKX = rmp.BKMK * x;
	//dm.BKMKur = rmp.BKMKJ * p_next + rmp.BKMKNX - rmp.BKMKX;
	//need these for volume preserving energy
	//Volume energy matrices How do I make thi
//	Eigen::MatrixXd BKMTIK = BT * fmp.KMTIK;
	//rmp.BKMTIH = BT * fmp.KMTIH ;
	//rmp.BKMTIKB = BKMTIK * B;
	//rmp.BKMTIG = BT* fmp.KMTIG;
	//rmp.BKMTIKJ = BKMTIK * J;
	//rmp.BKMTIKX = BKMTIK * x;
	

	rmp.BMJ = B.transpose() * fmp.M * J;

	{
		//Eigen::MatrixXd I = Eigen::MatrixXd::Identity(B.rows(), B.rows());
		
		Eigen::SparseMatrix<double> spI; spI.resize(B.rows(), B.rows());
		spI.setIdentity();
		
		Eigen::SparseMatrix<double> MJT, KMKJT, KMGT;
		MJT = fmp.MJ.transpose();
		KMKJT = fmp.KMKJ.transpose();
		KMGT = fmp.KMG.transpose();

		//rot_mat.IRI = LinearRotatedMatrix(I, I, W);
		rot_mat.B = LinearRotatedMatrix(spI, fmp.B_0, W);
		rot_mat.GmKB = LinearRotatedMatrix(fmp.GmK, B, W); //for arap/defo gradient
		rot_mat.BKMH = LinearRotatedMatrix(fmp.KMH.transpose(), B, W);
		rot_mat.BKMKJ = LinearRotatedMatrix(KMKJT, B, W);
		rot_mat.BKMKX = LinearRotatedMatrix(fmp.KMKX.transpose(), B, W);

	//	rot_mat.BKMKB = QuadraticRotatedMatrix(BT, B, fmp.KMK,  W); //quadratic term
		
	//	rot_mat.BMB = QuadraticRotatedMatrix(BT, B, fmp.M, W); //for inertia, quadratic term

	//	rot_mat.IRIRI = QuadraticRotatedMatrix(I, I, spI, W);
	
		rot_mat.BMJ = LinearRotatedMatrix(MJT, B, W); //for inertia
		rot_mat.BKMG = LinearRotatedMatrix(KMGT, B, W); //for volume grad.
	}

}

void FastCDSim::update_subspace(const Eigen::VectorXd& p, const Eigen::MatrixXd& W)
{
	Eigen::SparseMatrix<double> R;
	Eigen::MatrixXd R_stack;
	per_vertex_rotation(p, W, R, R_stack);

	B = R * fmp.B_0;

	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols()); //a few of these need rest space positions

	Eigen::MatrixXd BT = B.transpose();
	//rmp.BA = BT * fmp.A;
	//rmp.BAB = rmp.BA * B;
	//rmp.BAJ = rmp.BA * J;

	rmp.BM = BT * fmp.M;
	//rmp.BMB = rmp.BM * B;
	rmp.BMJ = rmp.BM * J;


	//need these for bending energy
	rmp.BKMK = BT * fmp.KMK;
	//rmp.BKMKB = rmp.BKMK * B;
	rmp.BKMH = BT * fmp.KMH;
	rmp.GmKB = fmp.GmK * B;
	rmp.BKMG = BT * fmp.KMG;
	rmp.BKMKJ = rmp.BKMK * J;
	rmp.BKMKX = rmp.BKMK * x;

}

void FastCDSim::update_subspace_fast(const Eigen::VectorXd& p)
{
	//All of these matrices require 
	Eigen::MatrixXd BKMH, BKMKX, GmKB, BKMKJ, BKMKB, BMB, BMJ, BKMG;


	//Eigen::SparseMatrix<double> R;
	//Eigen::MatrixXd R_stack;
//	per_vertex_rotation(p, W, R, R_stack);

	//Eigen::MatrixXd R_test;
	//rot_mat.IRI.update(p, R_test);
	//
	//Eigen::MatrixXd error = R_test - R;
	//printf("R test error : %f \n", error.cwiseAbs().sum());

	//B = R * fmp.B_0;

	rot_mat.B.update(p, B);
	//rmp.BMJ = B.transpose() * fmp.MJ;
    rot_mat.BMJ.update(p, BMJ); //for inertia
	//printf("BMJ error norm : %f \n", (BMJ.transpose() - rmp.BMJ).cwiseAbs().sum());
	rmp.BMJ = BMJ.transpose();
	
	//rmp.BKMH = B.transpose() * fmp.KMH;
	rot_mat.BKMH.update(p, BKMH);
    //printf("BKMH error norm : %f \n", (BKMH.row(0).transpose() - rmp.BKMH).cwiseAbs().sum());
	rmp.BKMH = BKMH.row(0).transpose();
	
	//rmp.GmKB =  fmp.GmK * B;
	rot_mat.GmKB.update(p, GmKB);
	//printf("GmKB error norm : %f \n",(GmKB - rmp.GmKB).cwiseAbs().sum());
	rmp.GmKB = GmKB;
	
	//
	//rmp.BKMG = B.transpose() * fmp.KMG;
	rot_mat.BKMG.update(p, BKMG); //for volume grad.
	//printf("BKMG error norm : %f \n", (BKMG.transpose() - rmp.BKMG).cwiseAbs().sum());
	rmp.BKMG = BKMG.transpose();
	//
	//rmp.BKMKJ = B.transpose() * fmp.KMKJ;
	rot_mat.BKMKJ.update(p, BKMKJ);
	//printf("BKMKJ error norm : %f \n", (BKMKJ.transpose() - rmp.BKMKJ).cwiseAbs().sum());
	rmp.BKMKJ = BKMKJ.transpose();
	//
	//rmp.BKMKX = B.transpose() * fmp.KMKX;
	rot_mat.BKMKX.update(p, BKMKX);
	//printf("BKMKX error norm : %f \n", (BKMKX.row(0).transpose() - rmp.BKMKX).cwiseAbs().sum());
	rmp.BKMKX = BKMKX.row(0).transpose();

	//rmp.BKMKB = B.transpose() * fmp.KMK * B;
	//rot_mat.BKMKB.update(p, BKMKB); //quadratic term
	//printf("BKMKB error norm : %f \n", (BKMKB.transpose() - rmp.BKMKB).cwiseAbs().sum());
	//rmp.BKMKB = BKMKB.transpose();

	//rmp.BMB = B.transpose() * fmp.M * B;
	//rot_mat.BMB.update(p, BMB); //for inertia, quadratic term
	//printf("BMB error norm : %f \n", (BMB - rmp.BMB).cwiseAbs().sum());
	//rmp.BMB = BMB.transpose();

//	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(B.rows(), B.rows());
//	Eigen::SparseMatrix<double> spI; spI.resize(B.rows(), B.rows());
//	spI.setIdentity();
//	Eigen::MatrixXd IRIRI = I * R.transpose() * I * R * I; Eigen::MatrixXd IRIRI_test;
//	rot_mat.IRIRI.update(p, IRIRI_test);
//	printf("IRIRI error norm : %f  \n", (IRIRI - IRIRI_test).norm());
}


void FastCDSim::init_clustering_matrices()
{
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	Eigen::SparseMatrix<double>  G_tmp, S_cols, S_rows;
	grouping_matrix_from_clusters(labels, fmp.G);

	//I wish there was a kron function
	G_tmp = igl::repdiag(fmp.G, 3);
	interweaving_matrix(fmp.G.cols(), 3, S_cols);
	interweaving_matrix(fmp.G.rows(), 3, S_rows);
	G_tmp = S_rows.transpose() * G_tmp * S_cols;
	G_tmp = igl::repdiag(G_tmp, 3);
	fmp.G_1 = G_tmp;

	Eigen::VectorXd cluster_mass;
	
	G_tmp = fmp.G * fmp.Vol;
	igl::sum(G_tmp, 2, cluster_mass);
	fmp.Vol_c = cluster_mass.asDiagonal();

	fmp.G_m = fmp.G_1 * fmp.Vol_exp;
	igl::sum(fmp.G_m, 2, cluster_mass);

	Eigen::SparseMatrix<double> cm, cmi;
	cm = cluster_mass.asDiagonal();
	igl::invert_diag(cm, cmi);
	fmp.G_m =  cmi* fmp.G_m  ;
	
	Eigen::VectorXd row_sum_check;
	igl::sum(fmp.G_m, 2, row_sum_check);
	//std::cout << row_sum_check << std::endl;

	//sm.G_exp * sm.FM * (sm.H + sm.K * u);
	//full space matrices, matrices needed for the full space quasi-newton gradient computation step
	
	{
		fmp.GmH = fmp.G_m * fmp.H;
		fmp.GmK = fmp.G_m * fmp.K;
		fmp.GmKJ = fmp.GmK * J;
		fmp.GmKX = fmp.GmK * x; 

		fmp.GMH = fmp.G_1* fmp.Vol_exp * fmp.H;
		fmp.GMK = fmp.G_1 * fmp.Vol_exp * fmp.K;
		fmp.KMG = fmp.K.transpose() * fmp.Vol_exp * fmp.G_1.transpose();//sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//	fmp.KMTIG = fmp.K.transpose() * fmp.Vol_exp * fmp.traceMat.transpose() * fmp.traceMat * fmp.G_1.transpose();

		fmp.GMKJ = fmp.GMK * J;
		fmp.GMKX = fmp.GMK * x;
	}
}



void FastCDSim::energy(const Eigen::VectorXd& u, const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev, double& bending, double& volume, double& inertia)
{
	const int l = do_clustering ? num_clusters : T.rows();
	Eigen::VectorXd FV_flat = fmp.GmH + fmp.GmK * u;//sm.G_exp * sm.FM * (sm.H + sm.K * u);
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::VectorXd volume_energy_tet = Eigen::VectorXd::Zero(l);

	Eigen::Matrix3d rot, F, cov, rot_cov;
	double tr;
	for (int i = 0; i < l; i++)
	{
		F = F_stack.block(3 * i, 0, 3, 3);
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;

		tr = (rot.transpose() * F).diagonal().sum() - 3.0;
		volume_energy_tet(i) = 0.5 * fmp.Vol_c.coeff(i, i) * tr * tr;
	}
	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

	Eigen::VectorXd u_next = u;
	Eigen::VectorXd y = (2.0 * u_curr) - u_prev;
	Eigen::VectorXd u_y = u - y;
	Eigen::VectorXd Mu_y = fmp.M * u_y;
	double uMu = u_y.transpose() * Mu_y;
	//std::cout << uMu << std::endl;
	//std::cout << 0.5 * (1.0 / (dt * dt));
	inertia = do_inertia * 0.5 * (1.0 / (dt * dt)) * uMu;

	Eigen::VectorXd f_r = fmp.H + fmp.K * u - fmp.G_1.transpose() * R_flat;

	bending = 0.5 * stiffness* f_r.transpose() * fmp.Vol_exp * f_r;
	//Eigen::VectorXd tr = fmp.traceMat * f_r;
	volume = 0.5 * incompressibility * (volume_energy_tet.sum());//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * (volume_energy_tet.sum());

	//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * volume;
	//printf("bending : %e, inertia : %e, volume : %e ", elastic, inertia, volume);
	//double e = elastic + inertia;
}


void FastCDSim::energy(const Eigen::VectorXd& uc_next, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev, const Eigen::VectorXd& ur_next, const Eigen::VectorXd& ur_curr, const Eigen::VectorXd& ur_prev, double& bending, double& volume, double& inertia)
{
	const int l = do_clustering ? num_clusters : T.rows();
	Eigen::VectorXd FV_flat = fmp.GmH + fmp.GmK * (uc_next + beta* ur_next);//sm.G_exp * sm.FM * (sm.H + sm.K * u);
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::VectorXd volume_energy_tet = Eigen::VectorXd::Zero(l);

	Eigen::Matrix3d rot, F, cov, rot_cov;
	double tr;
	for (int i = 0; i < l; i++)
	{
		F = F_stack.block(3 * i, 0, 3, 3);
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;

		tr = (rot.transpose() * F).diagonal().sum() - 3.0;
		volume_energy_tet(i) = 0.5 * fmp.Vol_c.coeff(i, i) * tr * tr;
	}
	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

	Eigen::VectorXd u_next = uc_next + ur_next;
	Eigen::VectorXd y = (2.0 * (uc_curr + ur_curr)) - (uc_prev + ur_prev);
	Eigen::VectorXd u_y = u_next - y;
	Eigen::VectorXd Mu_y = fmp.M * u_y;
	double uMu = u_y.transpose() * Mu_y;
	//std::cout << uMu << std::endl;
	//std::cout << 0.5 * (1.0 / (dt * dt));
	inertia = do_inertia * 0.5 * (1.0 / (dt * dt)) * uMu;

	Eigen::VectorXd f_r = fmp.H + fmp.K *(uc_next + beta*ur_next) - fmp.G_1.transpose() * R_flat;

	bending = 0.5 * stiffness * f_r.transpose() * fmp.Vol_exp * f_r;
	//Eigen::VectorXd tr = fmp.traceMat * f_r;
	volume = 0.5 * incompressibility * (volume_energy_tet.sum());//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * (volume_energy_tet.sum());

	//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * volume;
	//printf("bending : %e, inertia : %e, volume : %e ", elastic, inertia, volume);
	//double e = elastic + inertia;
}


void FastCDSim::energy_reduced(const Eigen::VectorXd& z, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, double& bending, double& volume, double& inertia)
{
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * 3 );
	Eigen::VectorXd u_curr, u_prev, u_next;
	Eigen::VectorXd ur_next = J * p_next - x;
	Eigen::VectorXd ur_curr = J * p_curr - x;
	Eigen::VectorXd ur_prev = J * p_prev - x;

	energy(B*z, B*z_curr, B*z_prev, ur_next, ur_curr, ur_prev, bending, volume, inertia);
}


void FastCDSim::energy_complementary_reduced(const Eigen::VectorXd& z, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, 
 double& bending, double& volume, double& inertia )
{

	Eigen::VectorXd u_next = B * z;
	Eigen::VectorXd u_curr = B * z_curr;
	Eigen::VectorXd u_prev = B * z_prev;

	energy(u_next, u_curr, u_prev, bending, volume, inertia);
}


void FastCDSim::kinetic_energy_complementary_full(const Eigen::VectorXd& uc, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev, double& inertia)
{
	Eigen::VectorXd u_next = uc;
	Eigen::VectorXd y = (2.0 * uc_curr) - uc_prev;
	Eigen::VectorXd u_y = uc - y;
	Eigen::VectorXd Mu_y = fmp.M * u_y;
	double uMu = u_y.transpose() * Mu_y;
	inertia =do_inertia *  0.5 * (1.0 / (dt * dt)) * uMu;
}

void FastCDSim::kinetic_energy_complementary_reduced(const Eigen::VectorXd& z, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev,double& inertia)
{

	Eigen::VectorXd uc = B * z, uc_curr = B * z_curr, uc_prev = B * z_prev;
	Eigen::VectorXd y = (2.0 * uc_curr) - uc_prev;
	Eigen::VectorXd u_y = uc - y;
	Eigen::VectorXd Mu_y = fmp.M * u_y;
	double uMu = u_y.transpose() * Mu_y;
	inertia = do_inertia *0.5 * (1.0 / (dt * dt)) * uMu;


}



void FastCDSim::grad(const Eigen::VectorXd& uc_next, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev, const Eigen::VectorXd& ur_next, const Eigen::VectorXd& ur_curr, const Eigen::VectorXd& ur_prev,
	Eigen::VectorXd& bending_grad, Eigen::VectorXd& volume_grad, Eigen::VectorXd& inertia_grad)
{
	const int l = do_clustering ? num_clusters : T.rows();

	Eigen::VectorXd FV_flat = (fmp.GmH + fmp.GmK * (uc_next + beta*ur_next));
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::Matrix3d rot, F;
	Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	//fit rotations
	for (int i = 0; i < l; i++)
	{
		F = (F_stack.block(3 * i, 0, 3, 3));
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;

		const double tr = (rot.transpose() * F).diagonal().sum() - 3;
		trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
	}
	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
	const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());

	Eigen::VectorXd y = (2.0 *( uc_curr + ur_curr) - (uc_prev- ur_prev));

	inertia_grad = do_inertia * (1.0 / (dt * dt)) * (fmp.M * ((uc_next + ur_next) - y));
	bending_grad = stiffness * fmp.KMH + fmp.KMK * (uc_next + beta*ur_next) - fmp.KMG * R_flat;
	volume_grad = incompressibility * fmp.KMG * gb_flat;

}
void FastCDSim::grad(const Eigen::VectorXd& u, const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev,  Eigen::VectorXd& bending_grad, Eigen::VectorXd& volume_grad, Eigen::VectorXd& inertia_grad)
{
	const int l = do_clustering ? num_clusters : T.rows();

	Eigen::VectorXd FV_flat = (fmp.GmH + fmp.GmK * u);
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::Matrix3d rot, F;
	Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	//fit rotations
	for (int i = 0; i < l; i++)
	{
		F = (F_stack.block(3 * i, 0, 3, 3));
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;

		const double tr = (rot.transpose() * F).diagonal().sum() - 3;
		trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
	}
	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
	const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());

	Eigen::VectorXd y = (2.0 * u_curr - u_prev);

	inertia_grad = do_inertia * (1.0 / (dt * dt)) * (fmp.M * ( u - y));
	bending_grad = stiffness * fmp.KMH + fmp.KMK * u - fmp.KMG * R_flat;
	volume_grad = incompressibility * fmp.KMG * gb_flat;
}

void FastCDSim::grad_reduced(const Eigen::VectorXd& z, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev,
 Eigen::VectorXd& bending_grad, Eigen::VectorXd& volume_grad, Eigen::VectorXd& inertia_grad)
{
	const int l = do_clustering ? num_clusters : T.rows();

	Eigen::VectorXd GmKur = fmp.GmKJ * p_next - fmp.GmKX;
	Eigen::VectorXd FV_flat = (fmp.GmH + rmp.GmKB * z + beta*GmKur);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::Matrix3d rot, F;


	//fit rotations
	for (int i = 0; i < l; i++)
	{
		F = (F_stack.block(3 * i, 0, 3, 3));
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;


		//calculate (tr R^T F -I)R
		const double tr = (rot.transpose() * F).diagonal().sum() - 3;
		trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
	}

	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
	const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
	Eigen::VectorXd BMy = rmp.BMB * (2.0 * z_curr - z_prev) + rmp.BMJ * (2.0 * p_curr - p_prev);  
	Eigen::VectorXd BMy_tilde = BMy - rmp.BMJ * p_next;

	inertia_grad = do_inertia * (1.0 / (dt * dt)) * (rmp.BMB * z - BMy_tilde);
	bending_grad = stiffness * (rmp.BKMH + rmp.BKMKB * z + beta* dm.BKMKur - rmp.BKMG * R_flat);
	volume_grad = incompressibility * rmp.BKMG * gb_flat;
}

void FastCDSim::grad_complementary(const Eigen::VectorXd& uc, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev,
	Eigen::VectorXd& bending_grad, Eigen::VectorXd& volume_grad, Eigen::VectorXd& inertia_grad)
{
	Eigen::VectorXd p = Eigen::VectorXd::Zero(J.cols()), p_curr = Eigen::VectorXd::Zero(J.cols()), p_prev = Eigen::VectorXd::Zero(J.cols());
	grad(uc, uc_curr, uc_prev,  bending_grad, volume_grad, inertia_grad);
}


void FastCDSim::grad_complementary_reduced(const Eigen::VectorXd& z, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, 
 Eigen::VectorXd& bending_grad, Eigen::VectorXd& volume_grad,	Eigen::VectorXd& inertia_grad)
{

	const int l = do_clustering ? num_clusters : T.rows();

	Eigen::VectorXd FV_flat = (fmp.GmH + rmp.GmKB * z);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
	Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
	Eigen::MatrixXd R_cov = R;

	Eigen::Matrix3d rot, F;


	//fit rotations
	for (int i = 0; i < l; i++)
	{
		F = (F_stack.block(3 * i, 0, 3, 3));
		igl::polar_svd3x3(F, rot);
		R.block(3 * i, 0, 3, 3) = rot;


		//calculate (tr R^T F -I)R
		const double tr = (rot.transpose() * F).diagonal().sum() - 3;
		trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
	}

	const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
	const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
	Eigen::VectorXd BMy = rmp.BMB * (2.0 * z_curr - z_prev) ;
	Eigen::VectorXd BMy_tilde = BMy;

	inertia_grad = do_inertia * (1.0 / (dt * dt)) * (rmp.BMB * z - BMy_tilde);
	bending_grad = stiffness * (rmp.BKMH + rmp.BKMKB * z + dm.BKMKur - rmp.BKMG * R_flat);
	volume_grad = incompressibility * rmp.BKMG * gb_flat;

}