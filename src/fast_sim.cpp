#include "fast_sim.h"
#include "deformation_gradient_from_u_prefactorized_matrices.h"
#include "trace_matrix_operator.h"
#include "grouping_matrix_from_clusters.h"
#include "interweaving_matrix.h"
#include "kmeans.h"
#include "compute_modes_matlab.h"
#include "compute_modes_spectra.h"
#include "compute_clusters_igl.h"
#include "compute_clusters_matlab.h"
#include "arap_hessian.h"

#include <igl/polar_svd3x3.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/volume.h>
#include <igl/sum.h>
#include <igl/colon.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/average_onto_faces.h>


#include <filesystem>

#include <igl/get_seconds.h>
//TODO: need to  implement all these methods... yikes
FastSim::FastSim(){}

FastSim::FastSim(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, const Eigen::SparseMatrix<double>& J, double ym, double pr, double dt, int num_modes, int num_clusters, std::string modes_file_dir, std::string clusters_file_dir, bool do_reduction, bool do_clustering, int num_modal_features)
: X(X), T(T), J(J), dt(dt), num_modes(num_modes), num_clusters(num_clusters), modes_file_dir(modes_file_dir), clusters_file_dir(clusters_file_dir), do_reduction(do_reduction), do_clustering(do_clustering), num_modal_features(num_modal_features)
{
	namespace fs = std::filesystem;


	stiffness = ym / (2.0 * (1.0 + pr));
	incompressibility = ym * pr / ((1.0 + pr) * (1.0 - 2.0 * pr));

	reduced_newton_solver = new DenseQuasiNewtonSolver();
	full_newton_solver = new QuasiNewtonSolver();

	do_clustering = false;
	do_reduction = false;

	if (!fs::exists(fs::path(modes_file_dir).parent_path()))
	{
		fs::create_directories(fs::path(modes_file_dir).parent_path());
	}
	if (!fs::exists(fs::path(clusters_file_dir).parent_path()))
	{
		fs::create_directories(fs::path(clusters_file_dir).parent_path());
	}


	init_modes(num_modes);
	init_clusters(num_clusters, num_modal_features);
	init_system_matrices();
	precompute_solvers();
}

Eigen::VectorXd FastSim::reduced_step(const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const  Eigen::VectorXd& bc) {
	
	if (!do_reduction)
	{
		std::cout << "simulator not configured for calling reduced_step(). Please switch simulator to reduced mode with sim.switch_reduced(true)..." << std::endl;
		return Eigen::VectorXd::Zero(0);
	}
	//dm.y = B * (2.0 * z_curr - z_prev);
	dm.BMy = rmp.BMB * (2.0 * z_curr - z_prev);            //u_curr and u_prev are total displacements, eg u_prev = uc_prev + ur_prev;
	dm.yMy = 4.0 * z_curr.transpose() * rmp.BMB * z_curr;
	dm.yMy += z_prev.transpose() * rmp.BMB * z_prev;
	dm.yMy -= 4.0 * z_curr.transpose() * rmp.BMB * z_prev;

	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> f;
	reduced_qnewton_energy_grad(f, grad_f);

	Eigen::VectorXd z_next;

	bool do_line_search = incompressibility < 1e-8 ? false : true;
	if (bc.rows() > 0)
	{
		z_next = reduced_newton_solver->solve_with_equality_constraints(z_curr, f, grad_f, bc, do_line_search);
	}
	else {
		z_next = reduced_newton_solver->solve(z_curr, f, grad_f, do_line_search);
	}
	return z_next;
}

Eigen::VectorXd FastSim::full_step(const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev,const  Eigen::VectorXd& bc)
{
	assert(J.rows() == bc.rows() && "Equality constraint matrix does not match rhs term");

	if (do_reduction)
	{
		std::cout << "simulator not configured for calling full_step. Please switch simulator to full mode with sim.switch_reduced(false)..." << std::endl;
		return Eigen::VectorXd::Zero(u_curr.rows());
	}
	dm.y = (2.0 * u_curr - u_prev);
	//dm.My = fmp.M * (2.0 * u_curr - u_prev);
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	full_qnewton_energy_grad(energy, grad_f);

	Eigen::VectorXd u_next;
	bool do_line_search = incompressibility < 1e-8 ? false : true;
	if (bc.rows() > 0) // should set a flag instead of doing this
	{
		u_next = full_newton_solver->solve_with_equality_constraints(u_curr, energy, grad_f, bc, do_line_search);
	}
	else
	{
		u_next = full_newton_solver->solve(u_curr, energy, grad_f, do_line_search);
	}

	return u_next;
}

Eigen::VectorXd FastSim::reduced_step_eval_metrics(const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& bc, double& e, double& grad_e)
{
	Eigen::VectorXd z_next = reduced_step(z_curr, z_prev, bc);
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	reduced_qnewton_energy_grad(energy, grad_f);

	e = energy(z_next);
	grad_e = grad_f(z_next).squaredNorm();
	return z_next;
}

Eigen::VectorXd FastSim::full_step_eval_metrics(const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev, const Eigen::VectorXd& bc, double& e, double& grad_e)
{
	Eigen::VectorXd u_next = full_step(u_curr, u_prev, bc);
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<double(const Eigen::VectorXd&)> energy;
	full_qnewton_energy_grad(energy, grad_f);

	e = energy(u_next); 
	grad_e = grad_f(u_next).squaredNorm();
	return u_next;
}

void FastSim::init_modes(int num_modes){

	namespace fs = std::filesystem;
	

	std::string B_file_path = modes_file_dir + "B.DMAT";
	std::string L_file_path = modes_file_dir + "L.DMAT";
	bool found_modes = igl::readDMAT(B_file_path, B_full);
	bool found_evals = igl::readDMAT(L_file_path, L_full);
	bool enough_modes = false;

	if (found_modes && found_evals)
	{
		printf("Found %i / %i cached modes at %s ... \n", B_full.cols(), num_modes, B_file_path.c_str());
		//std::cout << "Found cached modes at :" << modes_file_dir << std::endl;
		enough_modes = num_modes <= B_full.cols();

	}
	if (!enough_modes)
	{
	
		printf("Could not find cached modes at %s, computing %i modes from scratch...\n", B_file_path.c_str(), num_modes);
		//build constrained eigenvalue problem
		Eigen::SparseMatrix<double> H;
		Eigen::SparseMatrix<double> M;
		arap_hessian(X, T, H);
		H *= -1.0;
		//Eigen::SparseMatrix<double> I(M.rows(), M.cols());
		igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);

		//compute modes and load into S_full and B_full
		Eigen::MatrixXd modes;
#ifdef FAST_CD_USE_MATLAB
//#error "should not get here"
		compute_modes_matlab(H, M, num_modes, modes, L_full);
#else
		compute_modes_spectra(H, M, num_modes, modes, L_full);
#endif
		B_full = modes.block(0, 0, 3*X.rows(), num_modes);

		if (!fs::exists(fs::path(B_file_path).parent_path()))
		{
			fs::create_directories(fs::path(B_file_path).parent_path());
		}
		printf("Saving new modes at %s \n", B_file_path.c_str());
		//writeDMAt
		igl::writeDMAT(B_file_path, B_full, false);
		Eigen::MatrixXd L_full_mat = Eigen::Map<Eigen::MatrixXd>(L_full.data(), num_modes, 1); //annoying
		igl::writeDMAT(L_file_path, L_full_mat, false);
	}
	L = L_full.topRows(num_modes);
	B = B_full.block(0, 0, 3 * X.rows(), num_modes);
	//B.resize(X.rows() * 3, X.rows() * 3);
	//B.setIdentity();
}; //shared with parent

void FastSim::init_clusters(int num_clusters, int num_feature_modes)
{
	namespace fs = std::filesystem;
	
	const int l = do_clustering ? num_clusters : T.rows();
	std::string labels_file_path = clusters_file_dir + "labels_" + std::to_string(l) + "_features_" + std::to_string(num_feature_modes) + ".DMAT";

	bool found_clusters = igl::readDMAT(labels_file_path, labels);
	if (!found_clusters)
	{
		//10 modes is usually enough for a fine clustering
		if (! (l == T.rows()))
		{
			printf("Could not find cached clusters at %s,  clustering... \n", labels_file_path.c_str());
			double t_start = igl::get_seconds();
			Eigen::MatrixXd C;
			#ifdef FAST_CD_USE_MATLAB
						compute_clusters_matlab(T, B_full, L_full, num_clusters, num_feature_modes, labels, C);
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

}

void FastSim::init_system_matrices(){
	init_full_system_matrices();
	init_clustering_matrices();
	init_reduction_matrices();
} //overrides parent

void FastSim::init_full_system_matrices(){
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	Eigen::SparseMatrix<double> C;
	//cotan laplacian
	igl::cotmatrix(X, T, C);
	C = -igl::repdiag(C, 3);

	//mass matrix
	igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, fmp.M);
	fmp.M = igl::repdiag(fmp.M, 3);

	//constructs A matrix, and does preomputation if we will use it directly
	fmp.A = stiffness * C + (1.0 / (dt * dt))* fmp.M; //;

	//flattened_deformation_gradient matrices. 
	deformation_gradient_from_u_prefactorized_matrices(X, T, fmp.H, fmp.K, fmp.Vol_exp); 

	Eigen::VectorXd vol;
	igl::volume(X, T, vol);
	fmp.Vol = vol.asDiagonal();

	trace_matrix_operator(T.rows(), 3, fmp.traceMat);


	//Need these for bending energy
	fmp.KMH = fmp.K.transpose() * fmp.Vol_exp * fmp.H;
	fmp.KMK = fmp.K.transpose() * fmp.Vol_exp * fmp.K;				//this should be exactly equal to cotan laplacian

	fmp.HMH = fmp.H.transpose() * fmp.Vol_exp * fmp.H;



}

void FastSim::init_reduction_matrices()
{
	//init_modes(num_modes);
	rmp.BAB = B.transpose() * (fmp.A) * B;
	rmp.BMB = B.transpose() * (fmp.M) * B;
	rmp.BKMKB = B.transpose() * (fmp.KMK) * B;
	rmp.BKMH = B.transpose() * (fmp.KMH);
	rmp.BKMG = B.transpose() * (fmp.KMG);
	rmp.GmKB = fmp.GmK * B;

	rmp.SB = J * B;

}

void FastSim::init_clustering_matrices()
{
	//init_clusters(num_clusters, 10); //TODO expose this num_feature modes to UI
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	Eigen::SparseMatrix<double>  G, G_tmp, S_cols, S_rows;
	grouping_matrix_from_clusters(labels, G);

	//I wish there was a kron function
	G_tmp = igl::repdiag(G, 3);
	interweaving_matrix(G.cols(), 3, S_cols);
	interweaving_matrix(G.rows(), 3, S_rows);
	G_tmp = S_rows.transpose() * G_tmp * S_cols;
	G_tmp = igl::repdiag(G_tmp, 3);
	fmp.G_1 = G_tmp;

	Eigen::VectorXd cluster_mass;

	G_tmp = G * fmp.Vol;
	igl::sum(G_tmp, 2, cluster_mass);
	fmp.Vol_c = cluster_mass.asDiagonal();


	fmp.G_m = fmp.G_1 * fmp.Vol_exp;
	igl::sum(fmp.G_m, 2, cluster_mass);
	
	
	cluster_mass.array() = 1.0/cluster_mass.array();
	fmp.G_m = cluster_mass.asDiagonal() * fmp.G_m;

	
	//Eigen::VectorXd row_sum_check;
	//igl::sum(fmp.G_m, 2, row_sum_check);
	//std::cout << row_sum_check << std::endl;

	//sm.G_exp * sm.FM * (sm.H + sm.K * u);
	//full space matrices, matrices needed for the full space quasi-newton gradient computation step

	{

		fmp.GmH = fmp.G_m * fmp.H;
		fmp.GmK = fmp.G_m * fmp.K;

		fmp.GMH = fmp.G_1 * fmp.Vol_exp * fmp.H;
		fmp.GMK = fmp.G_1 * fmp.Vol_exp * fmp.K;

		fmp.KMG = fmp.K.transpose() * fmp.Vol_exp * fmp.G_1.transpose();//sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		

		fmp.GMG = fmp.G_1 * fmp.Vol_exp * fmp.G_1.transpose();
		fmp.HMG = fmp.H.transpose() * fmp.Vol_exp * fmp.G_1.transpose();


	}
}

void FastSim::full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
	f = [&](const Eigen::VectorXd u)
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
			volume_energy_tet(i) = 0.5*fmp.Vol_c.coeff(i, i) * tr * tr;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

		Eigen::VectorXd u_y = u - dm.y;
		double inertia = 0.5 * (1.0 / (dt * dt)) * u_y.transpose() * fmp.M * u_y;

		Eigen::VectorXd f_r = fmp.H + fmp.K * u - fmp.G_1.transpose() * R_flat;

		double bending = 0.5 * f_r.transpose() * fmp.Vol_exp * f_r;
		//Eigen::VectorXd tr = fmp.traceMat * f_r;
		//double volume = 0.5 * tr.transpose() * fmp.Vol * tr;//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		double elastic = stiffness * bending + incompressibility * (volume_energy_tet.sum());

		double e = elastic + inertia;
		return e;
	};

	grad_f = [&](const Eigen::VectorXd u)
	{
		const int l = do_clustering ? num_clusters : T.rows();
	
		Eigen::VectorXd FV_flat = (fmp.GmH + fmp.GmK * u); //get flattened, average deformation gradient
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());


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

		Eigen::VectorXd u_y = u - dm.y;
		//	Eigen::VectorXd traces = sm.traceMat *  sm.G_exp.transpose() * R_flat;
		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (fmp.M * u_y);
		Eigen::VectorXd bending_grad = fmp.KMH + fmp.KMK * u - fmp.KMG * R_flat;

		Eigen::VectorXd volume_grad = fmp.KMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad + incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + inertia_grad;
		return g;
	};

}

void FastSim::reduced_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
{
	f = [&](const Eigen::VectorXd z)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat = fmp.GmH + rmp.GmKB * z;//sm.G_exp * sm.FM * (sm.H + sm.K * u);
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd R_cov = R;

		Eigen::VectorXd volume_energy_tet = Eigen::VectorXd::Zero(l);

		Eigen::VectorXd bending_energy_tet = Eigen::VectorXd::Zero(l);
		Eigen::Matrix3d rot, F, cov, rot_cov;
		double tr, frob;
		for (int i = 0; i < l; i++)
		{
			F = F_stack.block(3 * i, 0, 3, 3);
			igl::polar_svd3x3(F, rot);
			R.block(3 * i, 0, 3, 3) = rot;

			tr = (rot.transpose() * F).diagonal().sum() - 3.0;
			volume_energy_tet(i) = 0.5 * fmp.Vol_c.coeff(i, i) * tr * tr;

			// too bad, assumes clustered F, but we don't do that for this term frob = ((F - R).transpose() * (F-R)).diagonal().sum();
			// too bad, assumes clustered F, but we don't do that for this term bending_energy_tet(i) = 0.5 * fmp.Vol_c.coeff(i, i)* frob * frob;
		}
		const Eigen::VectorXd R_flat = Eigen::Map<const Eigen::VectorXd>(R.data(), R.rows() * R.cols());

		double inertia = z.transpose() * rmp.BMB * z;
		inertia += -2.0 * z.transpose() * dm.BMy;
		inertia +=  dm.yMy;
		inertia *= 0.5 * (1.0 / (dt * dt));

		double bending = 0.5 * fmp.HMH + 0.5 * z.transpose() * rmp.BKMKB * z + z.transpose() * rmp.BKMH \
			 + 0.5 * R_flat.transpose() * fmp.GMG * R_flat \
		-fmp.HMG * R_flat -z.transpose() * rmp.BKMG * R_flat;


		//Eigen::VectorXd tr = fmp.traceMat * f_r;
		//double volume = 0.5 * tr.transpose() * fmp.Vol * tr;//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
		double elastic = stiffness * bending + incompressibility * volume_energy_tet.sum();

		double e = elastic + inertia;
		return e;
	};
	grad_f = [&](const Eigen::VectorXd z)
	{
		const int l = do_clustering ? num_clusters : T.rows();
		Eigen::VectorXd FV_flat = (fmp.GmH + rmp.GmKB * z);  //fmp.GMH + rmp.GMKB * z + dm.GMKur; //Perfect
		Eigen::MatrixXd F_stack = Eigen::Map<Eigen::MatrixXd>(FV_flat.data(), l * X.cols(), X.cols());
		Eigen::MatrixXd R = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());
		Eigen::MatrixXd trRT_IR = Eigen::MatrixXd::Zero(F_stack.rows(), F_stack.cols());

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
		Eigen::VectorXd inertia_grad = (1.0 / (dt * dt)) * (rmp.BMB*z - dm.BMy);
		Eigen::VectorXd bending_grad = rmp.BKMH + rmp.BKMKB*z - rmp.BKMG * R_flat;

		Eigen::VectorXd volume_grad = rmp.BKMG * gb_flat;
		Eigen::VectorXd elastic_grad = stiffness * bending_grad + incompressibility * volume_grad;
		Eigen::VectorXd g = elastic_grad + inertia_grad;

		return g;
	};
}

void FastSim::precompute_solvers() 
{
	reduced_newton_solver->precompute_with_equality_constraints(rmp.BAB, rmp.SB);
	//reduced_newton_solver->precompute_with_constraints(rmp.BTAB, constraint_prefactorization.SB);
	if (!do_reduction)
	{
			//Eigen::VectorXi all = igl::colon<int>(0, X.rows() * 3-1);
			//Eigen::VectorXi bi_flat = sm.J.cast<int>() * all;

		full_newton_solver->precompute_with_equality_constraints(fmp.A, J);
	}
}

//THESE ones are to update the sim in real time... dNot a priority
void FastSim::switch_clustering(bool do_clustering)
{
	this->do_clustering = do_clustering;
	init_clusters(num_clusters, num_modal_features);
	init_clustering_matrices();
	init_reduction_matrices();
}

void FastSim::switch_reduction(bool do_reduction) {

	this->do_reduction = do_reduction;
	init_modes(num_modes);
	init_system_matrices();
	precompute_solvers();
}

void FastSim::update_equality_constraint(const Eigen::SparseMatrix<double>& J)
{
	this->J = J;
	init_system_matrices();
	precompute_solvers();
}

void FastSim::update_material_properties(double ym, double pr)
{
	stiffness = ym / (2.0 * (1.0 + pr));
	incompressibility = ym * pr / ((1.0 + pr) * (1.0 - 2.0 * pr));

	init_system_matrices();

	precompute_solvers();
}

void FastSim::update_timestep(double dt)
{
	this->dt = dt;
	init_system_matrices();
	precompute_solvers();
}

void FastSim::update_modes(double new_num_modes)
{
	num_modes = new_num_modes;

	//recompute modes if we don't already have them loaded in (don't bother looking in cache)
	if (new_num_modes <= B_full.cols())
	{
		L = L_full.topRows(num_modes);
		B = B_full.block(0, 0, 3*X.rows(), num_modes);
	}
	else
	{
		init_modes(num_modes);
	}

	init_system_matrices();


	precompute_solvers();
}

void FastSim::update_clusters(double new_num_clusters)
{
	num_clusters = new_num_clusters;
	init_clusters(num_clusters, num_modal_features);
	init_system_matrices();

}

void FastSim::update_modes_cache_dir(std::string new_mode_dir)
{
	modes_file_dir = new_mode_dir;
}

void FastSim::update_clusters_cache_dir(std::string new_clusters_dir)
{
	clusters_file_dir = new_clusters_dir;
}


void FastSim::energy(const Eigen::VectorXd& u, double& bending, double& volume, double& inertia)
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

	Eigen::VectorXd u_y = u - dm.y;
	inertia = 0.5 * (1.0 / (dt * dt)) * u_y.transpose() * fmp.M * u_y;

	Eigen::VectorXd f_r = fmp.H + fmp.K * u - fmp.G_1.transpose() * R_flat;

	bending = 0.5 * f_r.transpose() * fmp.Vol_exp * f_r;
	//Eigen::VectorXd tr = fmp.traceMat * f_r;
	volume = 0.5 * incompressibility * (volume_energy_tet.sum());//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * (volume_energy_tet.sum());

	//sm.K.transpose()* sm.FM* sm.traceMat.transpose()* sm.traceMat* (sm.H + sm.K * u - sm.G_exp.transpose() * R_flat); //sm.K.transpose() * sm.FM * (sm.H + sm.K * u - sm.G_exp.transpose() * (R_flat)) ; //every tet in each cluster has the same
	//elastic = stiffness * bending + incompressibility * volume;
	//printf("bending : %e, inertia : %e, volume : %e ", elastic, inertia, volume);
	//double e = elastic + inertia;
}