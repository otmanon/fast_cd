#pragma once
#include <Eigen\Core>
#include <Eigen\Sparse>
#include "FastCDLocalGlobalSolver.h"
#include "CDLocalGlobalSolver.h"
#include "CDNewtonSolver.h"
#include "FastCDNewtonSolver.h"
#include "fast_sim.h"
/// <summary>
/// This class implements our core fast Complementary dynamics algorithm. It sets up the reduced and clustered physical system
/// and allows us to query a sample next step complementary displacement, taking reduced parameters as input.
/// 
/// This class by design does not store the simulation state. This is left up to the caller to give us the simulation state as desired.
/// </summary>
class FastCDSim : public FastSim
{
public:
	FastCDSim() {};

	/*
	Constructor for a FastCDSimulator that we can use to step complementary dynamics simulations, with our reduced method AND our the full method
	Input:
	X -	|V| x 3 matrix of mesh vertices
	T - |T|x4   tet mesh indices
	J - 3V x p  rig jacobian, can easily be obtained from the Rig class
	ym - double young's modulus, physical material parameter denoting stiffness
	pr - double poisson ratio, ranges from 0-0.49  (anything over 0.5, will explode). Controlls material incompressibility
	dt - double timestep
	num_modes - int  the number of modes to use for our accelerated CD scheme
	num_clusters - int the number of clusters to use for our accelerated CD scheme
	modes_file_dir - std::string - file directory where to find cached B.DMAT file AND the S.DMAT file for our modes, denoting eigenvectors and corresponding eigenvalues respectively. (Will recompute if unfindable)
	clusters_file_dir - std::string - file directory where to find cached labels_<num_labels>_features_<num_features_used_for_labels>.DMAT file, containing a label for each tet denoting which cluster it belongs to (Will recompute if unfindable)
	modal_features (optional) int - the number of modal features to use in our clustering scheme. We claim that number of features should have diminishing returns.
	*/
	FastCDSim(Eigen::MatrixXd& X,  Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& J, double ym, double pr, double dt, int num_modes, int num_clusters,  std::string modes_file_dir, std::string clusters_file_dir, bool do_reduction=true, bool do_clustering=true, int num_modal_features=10 );

	void init_modes(int num_modes);
	
	void init_clusters(int num_clusters, int num_feature_modes);
	
	void init_full_system_matrices();
	
	void init_reduction_matrices();

	void init_clustering_matrices();


	void full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

	void reduced_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

	/*
	Builds the constraint matrix C, and the corresponding rhs vector f and adds them to the system we solve
	Inputs: bI - c x 1 vertex indices in our mesh
			bc - c x 3 desired vertex positions

	*/
	void make_positional_constraints(Eigen::VectorXi& bI, Eigen::MatrixXd& bc);

	void update_positional_constraints(Eigen::MatrixXd& bc);

	Eigen::VectorXd reduced_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev);

	Eigen::VectorXd full_step(const Eigen::VectorXd& p_next, const Eigen::VectorXd& p_curr, const Eigen::VectorXd& p_prev, const Eigen::VectorXd& uc_curr, const Eigen::VectorXd& uc_prev);


	void update_compelementary_constraint(const Eigen::SparseMatrix<double>& J, std::string new_modes_dir, std::string new_clusters_dir);


	void precompute_solvers();

public:
	//Rig jacobian
	Eigen::SparseMatrix<double> J;

private:

	

	struct ConstraintMatrices{

		//Sparse equality matrix for positional constraints
		Eigen::SparseMatrix<double> S;
		//Desired user constraints used in the system solve. Deformed positional boundary conditions.
		Eigen::VectorXd Pbc;

		//same as the  before, but for reduced space constraints
		Eigen::MatrixXd SB;
		Eigen::MatrixXd SJ, SX;

	}constraint_prefactorization;



	struct FullSimPrecomputedMatrices
	{
		//full space hessian approximation (for ARAP: stiffness * L - (1/dt)^2 M)
		Eigen::SparseMatrix<double> A;

		//full space cotan Laplacian TODO: get rid of this, we may not need it anymore.
		Eigen::SparseMatrix<double> C;

		//complementary equality constraint J^T M D from complementary dynamics paper
		Eigen::SparseMatrix<double> Aeq;

		// diagonal mass matrix |V|x|V|
		Eigen::SparseMatrix<double> M;

		//grouping matrix that sums a scalar value across clusters |#clusters|x |T|, followed by the same thing but expanded out for all 9 entries in R
		Eigen::SparseMatrix<double> G, G_exp, G_m;

		//precomputed constant multiplications, needed for inertia, stuff like that
		Eigen::SparseMatrix<double> MJ;
		Eigen::VectorXd MX;

		// dF/du matrix. Can be derived by relating F to u: vec(F) = Kvec(u) + H
		Eigen::SparseMatrix<double> K;
		Eigen::VectorXd H; // constant term in using K to map from complementary displacements to deformation gradients ie F = H + K*u
		Eigen::SparseMatrix<double> traceMat; //trace linear operator
		Eigen::SparseMatrix<double> Vol; //volume for each tet.
		Eigen::SparseMatrix<double> Vol_exp; //expanded volume for each tet that multiplies each entry in flattened vec(F) by the tet's mass.
	
		//Bending energy (ARAP) matrices
		//matrices used in full space simulation (included grouping matrix by default
		Eigen::SparseMatrix<double> KMK, KMG, GMK;
		//Maps from full tets to clustered tets
		Eigen::SparseMatrix<double> GK, GFMK;

		Eigen::SparseMatrix<double> KMKJ;
		Eigen::VectorXd GMH, KMH, KMKNX, KMKX;
		Eigen::VectorXd GmH, GmKX;
		Eigen::SparseMatrix<double> GmK, GmKJ;
		Eigen::MatrixXd GMKJ;
		Eigen::VectorXd GMKNX, GMKX;

		//Volume energy matrices
		//precomputed matrices for incompresisble part of full simulation. Has a trace term in the middle
		Eigen::VectorXd KMTIH;
		Eigen::SparseMatrix<double> KMTIK;
		Eigen::SparseMatrix<double> KMTIG;


	
	}fmp;

	struct ReducedSimPrecomputedMatrices
	{
		//reduced order matrices
		Eigen::MatrixXd BTAB, BTA, BTAJ, BTMB, BTM, BTMJ;
		Eigen::VectorXd BTMX;//use this incase of BTMJ in the event of a null rig

		//Bending energy matrices
		//same as above but for reduced space simulation
		Eigen::MatrixXd GMKB, KB;
		Eigen::MatrixXd GmKB;
		Eigen::MatrixXd BMB;
		Eigen::MatrixXd BKMK, BKMKB, BKMG, BKMGm;
		Eigen::VectorXd BKMH;

		Eigen::MatrixXd BKMKJ;
		Eigen::VectorXd BKMKNX, BKMKX;

		//Volume energy matrices
		Eigen::MatrixXd  BKMTIKB, BKMTIG, BKMTIKJ;
		Eigen::VectorXd BKMTIH, BKMTIKX, BKMTIKNX;


	}rmp;
	
	

	struct dynamic_matrices
	{

		Eigen::VectorXd r;
		Eigen::VectorXd y_tilde;
		Eigen::VectorXd My_tilde;
		Eigen::VectorXd BMy_tilde;
		Eigen::VectorXd Kur;
		Eigen::VectorXd GMKur;
		Eigen::VectorXd KMKur;
		Eigen::VectorXd BKMKur;
		Eigen::VectorXd BKMTIKur;

		Eigen::VectorXd GmKur;
	} dm;


	CDNewtonSolver* full_newton_solver;
	FastCDNewtonSolver* reduced_newton_solver;


};