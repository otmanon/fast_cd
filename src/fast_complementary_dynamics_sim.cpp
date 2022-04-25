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
#include "kmeans.h"
#include "arap_hessian.h"
#include "trace_matrix_operator.h"
#include "covariance_scatter_matrix.h"
#include "deformation_gradient.h"
#include "interweaving_matrix.h"
#include "grouping_matrix_from_clusters.h"
#include "fast_complementary_dynamics_constrained_gevp.h"
#include "complementary_equality_constraint.h"
#include "compute_modes_matlab.h"
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "igl/sum.h"
#include "igl/count.h"
#include <igl/colon.h>
#include <igl/volume.h>
#include <igl/get_seconds.h>

#ifdef WIN32
	#include <filesystem>
#else
	#include <experimental/filesystem>
#endif

FastCDSim::FastCDSim( Eigen::MatrixXd& X,  Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& J, double ym, double pr, double dt, int num_modes, int num_clusters, std::string modes_file_dir, std::string clusters_file_dir, bool do_reduction, bool do_clustering, int num_modal_features)
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

	this->num_modal_features = num_modal_features;
	stiffness = ym / (2.0 * (1.0 + pr));
	incompressibility = ym* pr / ((1.0 + pr) * (1.0 - 2.0 * pr));



	Eigen::VectorXi _n; //throwaway var
	igl::boundary_facets(T, F, FiT, _n);
	//get list of exterior vertex indices
	igl::unique(F, ext_ind);
	
	reduced_newton_solver = new FastCDNewtonSolver(10, 1e-10);
	full_newton_solver = new CDNewtonSolver(10, 1e-6, 1, 1);

	update_compelementary_constraint(J, modes_file_dir, clusters_file_dir);
}


Eigen::VectorXd FastCDSim::reduced_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev)
{
	if (!do_reduction)
	{
		std::cout << "simulator not configured for calling reduced_step(). Please switch simulator to reduced mode with sim.switch_reduced(true)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	Eigen::VectorXd r, ur, BTMy;
	//have a valid rig

	BTMy = rmp.BTMB *  (2.0 * z_curr - z_prev) + rmp.BTMJ * (2.0*p_curr - p_prev);            //u_curr and u_prev are total displacements, eg u_prev = uc_prev + ur_prev;

	dm.BMy_tilde = BTMy - rmp.BTMJ * p_next;


	//dm.Kur = fmp.K * J * p_next - fmp.K * x;	
//	dm.GMKur = fmp.G_exp* fmp.Vol_exp* fmp.K* ( J* p_next - x);
	dm.GMKur = fmp.GMKJ * p_next   -fmp.GMKX;
//	dm.KMKur = fmp.KMKJ * p_next + fmp.KMKNX - fmp.KMKX;
	dm.BKMKur = rmp.BKMKJ * p_next  - rmp.BKMKX;
	
	dm.BKMTIKur = rmp.BKMTIKJ * p_next - rmp.BKMTIKX;
	dm.GmKur = fmp.GmKJ * p_next  - fmp.GmKX;

	Eigen::VectorXd uc_bc;

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> f;
	reduced_qnewton_energy_grad(f, grad_f);

	Eigen::VectorXd z_next;
	
	z_next = reduced_newton_solver->solve(z_curr, f, grad_f);
	
	return z_next;

}

Eigen::VectorXd FastCDSim::full_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev)
{
	if (do_reduction)
	{
		std::cout << "simulator not configured for calling full_step. Please switch simulator to full mode with sim.switch_reduced(false)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	Eigen::VectorXd My, r;
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	r = J * p_next + N*x;
	//compute dynamic matrix at the start of the timestep
	dm.r = r;//TODO: delete this once we're done playing around with it... unnecessary
	//precompute_with_constraints AJ and MJ for more speed
	My = fmp.M * (2.0 *  uc_curr -  uc_prev) + fmp.MJ * (2.0 *  p_curr - p_prev) + fmp.M *N * x;
	

	dm.My_tilde = My - fmp.M * dm.r;// +2 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	dm.GmKur = fmp.GmKJ * p_next + fmp.GmK * N * x - fmp.GmKX;
//	dm.y_tilde = y - dm.r;// +2.0 * sm.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	full_qnewton_energy_grad(energy, grad_f);

	Eigen::VectorXd uc_bc;
	Eigen::VectorXd uc_next;
	
	uc_next = full_newton_solver->solve(uc_curr, energy,  grad_f);
	


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
	init_system_matrices();
	
	precompute_solvers();
}

void FastCDSim::precompute_solvers()
{
	
	reduced_newton_solver->precompute(rmp.BTAB);
	if (!do_reduction)
	{
		full_newton_solver->precompute(fmp.A, fmp.Aeq);
	}
	
}

void FastCDSim::reduced_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
	f = [&](const Eigen::VectorXd z)
	{
		/*
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat = sm.GMH + sm.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F, cov, rot_cov;
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		double inertia = 0.5*(1.0 / (dt * dt)) * z.transpose() * (sm.BMB * z - 2.0*dm.BMy_tilde);

		Eigen::VectorXd F_R = (sm.H + sm.K * (B * z + dm.r) ) - sm.G_exp.transpose() * R_flat;
		double bending = 0.5*F_R.transpose() * sm.FM * F_R;
		double e = stiffness * bending + inertia;

		printf("bending : %e, inertia : %e \n ", bending, inertia);
		return e;*/

		/*
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd uc = B * z;
		//Eigen::VectorXd u = uc + dm.r - Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
		Eigen::VectorXd FV_flat = sm.GMH + sm.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F, cov, rot_cov;
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		//	Eigen::VectorXsd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;

		double inertia = 0.5 * (1.0 / (dt * dt)) * z.transpose() * (sm.BMB * z - 2.0 * dm.BMy_tilde);
		//	inertia +=  0.5* (1.0 / (dt * dt)) *dm.y_tilde.transpose() * sm.M * dm.y_tilde;
		Eigen::VectorXd f_r = sm.H + sm.KB * z + dm.Kur - sm.G_exp.transpose() * R_flat;
		//double bending = 0.5* f_r.transpose();
		double bending = 0.5 * f_r.transpose() * sm.FM * f_r;
		//Eigen::VectorXd tr = sm.traceMat * f_r;
		//double volume = 0.5 * tr.transpose() * sm.Vol * tr;//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		Eigen::VectorXd a = sm.H + dm.Kur;
		Eigen::VectorXd aC = a.transpose() * (sm.G_exp.transpose())
		//	double elastic = stiffness * bending;// +0 * incompressibility * volume;
		printf("bending : %e, inertia : %e \n ", elastic, inertia);
		double e = elastic + inertia;
		return e;*/
		return 0;

	};
	grad_f = [&](const Eigen::VectorXd z)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat =  (fmp.GmH + rmp.GmKB * z + dm.GmKur);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F;
		//fit rotations
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;

			//calculate (tr R^T F -I)R
			const double tr = (rot.transpose() * F).diagonal().sum() - 3;
			trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
		//	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;
		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB * z - dm.BMy_tilde);
		//Eigen::VectorXd bending_grad =  rmp.BKMG * (FV_flat - R_flat);
		
		Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB * z + dm.BKMKur - rmp.BKMG * R_flat;
		//Eigen::VectorXd volume_grad = rmp.BKMTIH +  rmp.BKMTIKB* z + dm.BKMTIKur  - rmp.BKMTIG * R_flat; //sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		Eigen::VectorXd volume_grad = rmp.BKMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad +incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + inertia_grad;
		return g;
		/*
		const int l = do_clustering ? num_clusters : T.rows();

		Eigen::VectorXd FV_flat = fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F, cov, rot_cov;
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3).transpose();
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot.transpose();
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB * z - dm.BMy_tilde);
		Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB *z + dm.BKMKur  - rmp.BKMG * R_flat;//sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
																							//by the way K^T M K is the same as the cotan laplacian!!
		Eigen::VectorXd volume_grad = rmp.BKMTIH+ rmp.BKMTIKB * z + dm.BKMTIKur - rmp.BKMTIG * R_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad + 0*incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + inertia_grad;

		return g;*/
	};
}

void FastCDSim::full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{

	f = [&](const Eigen::VectorXd uc)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd u = uc + dm.r - Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

		Eigen::VectorXd FV_flat = fmp.GMH + fmp.GMK * u;//sm.G_exp * sm.FM * (sm.H + sm.K * u);
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F, cov, rot_cov;
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		//	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;

		double inertia = 0.5 * (1.0 / (dt * dt)) * uc.transpose() * (fmp.M * uc - 2.0*dm.My_tilde);

		Eigen::VectorXd f_r = fmp.H + fmp.K * u - fmp.G_exp.transpose() * R_flat;

		//double bending = 0.5* f_r.transpose();
		double bending = 0.5* f_r.transpose() * fmp.Vol_exp * f_r;
		Eigen::VectorXd tr = fmp.traceMat * f_r;
		double volume = 0.5* tr.transpose() * fmp.Vol * tr;//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		double elastic = stiffness * bending +  incompressibility * volume;
		//printf("bending : %e, inertia : %e, volume : %e ", elastic, inertia, volume);
		double e = elastic + inertia;
		return e;
	};

	grad_f = [&](const Eigen::VectorXd uc)
	{
		const int l = do_clustering ? num_clusters : T.rows();
        Eigen::VectorXd u = uc + dm.r - Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

		Eigen::VectorXd FV_flat = (fmp.GmH + fmp.GmK * u);//fmp.GMH + fmp.GMK*u;//sm.G_exp * sm.FM * (sm.H + sm.K * u);
        Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::Matrix3d rot, F; 
		Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		//fit rotations
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;

			const double tr = (rot.transpose() * F).diagonal().sum() - 3;
			trRT_IR.block(3 * i, 0, 3, 3) = tr * rot;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());
		const Eigen::VectorXd gb_flat = Eigen::Map<const Eigen::VectorXd>(trRT_IR.data(), R.rows() * R.cols());
	//	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;
		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (fmp.M *uc -  dm.My_tilde);
		Eigen::VectorXd bending_grad = fmp.KMH + fmp.KMK * u - fmp.KMG * R_flat;
	//	Eigen::VectorXd bending_grad =  fmp.KMG * (FV_flat - R_flat);
	//	Eigen::VectorXd volume_grad = fmp.KMTIH + fmp.KMTIK * u - fmp.KMTIG * R_flat; //sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		Eigen::VectorXd volume_grad = fmp.KMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad +incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + inertia_grad;
		return g;
	};

}

/*
Makes sure the simulation enforces positional constraints on mesh indices bi. bc is a list of desired vertex positions which will pe 
enforce as hard constraints in the optimizaiton
*/
void FastCDSim::make_positional_constraints(Eigen::VectorXi& bi, Eigen::MatrixXd& bc)
{
	
	//Need to make constraint index matrix C  ! This is 3n x 3c matrix This is essentially a selection matrix that maps each constraint to its vertex
	std::vector<Eigen::Triplet<double >> tripletList;
	tripletList.reserve(bi.rows() * 3);
	for (int i = 0; i < bi.rows(); i++)
	{	
		tripletList.emplace_back(Eigen::Triplet<double>( i + bi.rows() * 0, bi(i) + X.rows() * 0, 1));
		tripletList.emplace_back(Eigen::Triplet<double>(i + bi.rows() * 1, bi(i) + X.rows() * 1, 1));
		tripletList.emplace_back(Eigen::Triplet<double>(i + bi.rows() * 2, bi(i) + X.rows() * 2, 1));
	}
	constraint_prefactorization.S.resize( 3*bi.rows(), 3 * X.rows());
	constraint_prefactorization.S.setFromTriplets(tripletList.begin(), tripletList.end());

	constraint_prefactorization.Pbc = Eigen::Map<Eigen::VectorXd>(bc.data(), bc.rows() * bc.cols()); // additionally will have to change this each timestep because we are solving for displacement, not position
																// this is constrained position = uc + ur + x. We solve_with_constraints for uc so constraint should become Cu^c = f - Cr													
	//precompute some constraint matrices.
	constraint_prefactorization.SB = constraint_prefactorization.S * B;
	constraint_prefactorization.SJ = constraint_prefactorization.S * J;
	constraint_prefactorization.SX = constraint_prefactorization.S * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows()*X.cols());

	precompute_solvers();
}


void FastCDSim::update_positional_constraints(Eigen::MatrixXd& bc)
{
	assert(constraint_prefactorization.S.rows() > 0 && "Cannot update physical constraints if we did not make them first. Please provide list of constrained indices \
								and call FastCDSim.make_positional_constraints(indeces, constraints)");
	constraint_prefactorization.Pbc = Eigen::Map<Eigen::VectorXd>(bc.data(), bc.rows() * bc.cols());
}



void FastCDSim::init_modes(int num_modes)
{

#ifdef WIN32
		namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif
	std::string B_file_path = modes_file_dir + "B.DMAT";
	std::string L_file_path = modes_file_dir + "S.DMAT";
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
		fast_complementary_dynamics_constrained_geivp(X, T, J, H, M);

		//compute modes and load into S_full and B_full
		Eigen::MatrixXd modes;
		compute_modes_matlab(H, M, num_modes, modes, L_full);
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
}

void FastCDSim::init_clusters(int num_clusters, int num_feature_modes)
{

	#ifdef WIN32
		namespace fs = std::filesystem;
	#else
		namespace fs = std::experimental::filesystem;
	#endif
	const int l = do_clustering ? num_clusters : T.rows();
	std::string labels_file_path = clusters_file_dir + "labels_" + std::to_string(l) + "_features_" + std::to_string(num_feature_modes) + ".DMAT";

	bool found_clusters = igl::readDMAT(labels_file_path, labels);
	if (!found_clusters)
	{
		if (!(l == T.rows()))
		{
			printf("Could not find cached clusters at %s,  clustering...\n", labels_file_path.c_str());
			//10 modes is usually enough for a fine clustering
			num_feature_modes = num_feature_modes < B_full.cols() ? num_feature_modes : B_full.cols();
			Eigen::RowVectorXd L = L_full.topRows(num_feature_modes).transpose();
			Eigen::MatrixXd B = B_full.block(0, 0, 3 * X.rows(), num_feature_modes);
			Eigen::MatrixXd B_verts = Eigen::Map<Eigen::MatrixXd>(B.data(), B.rows() / 3, B.cols() * 3);
			//optionally normalize by S
			double t_start = igl::get_seconds();
			B.array().rowwise() /= L.row(0).array();
			B.array().rowwise() /= L.row(0).array();
			Eigen::MatrixXd B_faces, C;
			igl::average_onto_faces(T, B_verts, B_faces);
			igl::kmeans(B_faces, num_clusters, C, labels);
			printf("Done clustering! took %g seconds \n", igl::get_seconds() - t_start);

			if (!fs::exists(fs::path(labels_file_path).parent_path()))
			{
				fs::create_directories(fs::path(labels_file_path).parent_path());
			}
			printf("Saving at %s...\n", labels_file_path.c_str());
			Eigen::MatrixXi labels_mat = Eigen::Map<Eigen::MatrixXi>(labels.data(), labels.rows(), 1);
			igl::writeDMAT(labels_file_path, labels_mat, false);
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
	fmp.A = stiffness * fmp.C + (1.0 / dt) * (1.0 / dt) * fmp.M;

	complementary_equality_constraint(X, T, J, fmp.Aeq);
	fmp.MJ = fmp.M * J;
	fmp.MX = fmp.M * Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());

	//flattened_deformation_gradient matrices. 
	deformation_gradient_from_u_prefactorized_matrices(X, T, fmp.H, fmp.K, fmp.Vol_exp);

	Eigen::VectorXd vol;
	igl::volume(X, T, vol);
	fmp.Vol = vol.asDiagonal();
	trace_matrix_operator(T.rows(), 3, fmp.traceMat);
	
	
	//Need these for bending energy
	fmp.KMH = fmp.K.transpose() * fmp.Vol_exp * fmp.H;
	fmp.KMK = fmp.K.transpose() * fmp.Vol_exp * fmp.K;				//this should be exactly equal to cotan laplacian
	
	fmp.KMTIH = fmp.K.transpose() * fmp.Vol_exp * fmp.traceMat.transpose() * fmp.traceMat * fmp.H;
	fmp.KMTIK = fmp.K.transpose() * fmp.Vol_exp * fmp.traceMat.transpose() * fmp.traceMat * fmp.K;

	//should do clustering matrices here.

	fmp.KMKJ = fmp.KMK * J; 
	fmp.KMKX = fmp.KMK * x;



	//build H
	// build K
}

void FastCDSim::init_reduction_matrices()
{
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols()); //a few of these need rest space positions
	//general system matrices
	rmp.BTA = B.transpose() * fmp.A;
	rmp.BTAB = rmp.BTA * B;
	rmp.BTAJ = rmp.BTA * J;

	rmp.BTM = B.transpose() * fmp.M;
	rmp.BTMB = rmp.BTM * B;
	rmp.BTMJ = rmp.BTM * J;

	rmp.BTMX = rmp.BTM * x;

	//need these for bending energy
	rmp.BKMK = B.transpose() * fmp.KMK;
	rmp.BKMKB = rmp.BKMK * B;
	rmp.BMB = rmp.BTMB; //same thing 
	rmp.BKMH = B.transpose() * fmp.KMH;

	rmp.GMKB = fmp.GMK * B;
	rmp.GmKB = fmp.GmK * B;
	rmp.BKMG = B.transpose() * fmp.KMG;

	rmp.BKMGm = B.transpose() * (fmp.K.transpose() * (fmp.Vol_exp * fmp.G_m.transpose()));
	rmp.BKMKJ = rmp.BKMK * J;
	rmp.BKMKX = rmp.BKMK * x;
	//dm.BKMKur = rmp.BKMKJ * p_next + rmp.BKMKNX - rmp.BKMKX;
	//need these for volume preserving energy
	//Volume energy matrices How do I make thi
	Eigen::MatrixXd BKMTIK = B.transpose() * fmp.KMTIK;
	rmp.BKMTIH = B.transpose() * fmp.KMTIH ;
	rmp.BKMTIKB = BKMTIK * B;
	rmp.BKMTIG = B.transpose() * fmp.KMTIG;
	rmp.BKMTIKJ = BKMTIK * J;
	rmp.BKMTIKX = BKMTIK * x;
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
	fmp.G_exp = G_tmp;

	Eigen::VectorXd cluster_mass;
	
	fmp.G_m = fmp.G_exp * fmp.Vol_exp;
	igl::sum(fmp.G_m, 2, cluster_mass);
	fmp.G_m =  cluster_mass.asDiagonal().inverse()* fmp.G_m  ;
	
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

		fmp.GMH = fmp.G_exp * fmp.Vol_exp * fmp.H;
		fmp.GMK = fmp.G_exp * fmp.Vol_exp * fmp.K;
		fmp.KMG = fmp.K.transpose() * fmp.Vol_exp * fmp.G_exp.transpose();//sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		fmp.KMTIG = fmp.K.transpose() * fmp.Vol_exp * fmp.traceMat.transpose() * fmp.traceMat * fmp.G_exp.transpose();


		fmp.GMKJ = fmp.GMK * J;
		fmp.GMKX = fmp.GMK * x;
	}
}

