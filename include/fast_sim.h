#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "QuasiNewtonSolver.h"
#include "DenseQuasiNewtonSolver.h"
//TODO: Make THIS the parent class of fast simulators, and have FastCD Sim inherit from it
class FastSim
{
public:
	//TODO: need to  implement all these method
	FastSim();

	/*
		Creates a simulator that is controlled by linear equality constraints S. (In our case get this matrix from a rig)
		X -  nx3 rest mesh pose
		T -  |T|x4 tet indices into X
		S -  (c x 3n) linear equality constraints on flattened X such that Su = f specifies vertex displacements of f for all constrained vertices

		ym - youngs modulus
		pr - poisson ratio
		dt - timestep
		*/
	 FastSim(const Eigen::MatrixXd& X,const  Eigen::MatrixXi& T, const Eigen::SparseMatrix<double>& S, double ym, double pr, double dt, int num_modes, int num_clusters, std::string modes_file_dir, std::string clusters_file_dir, bool do_reduction = true, bool do_clustering = true, int num_modal_features=10);
	
	 /*
	 * NOT A PRIORITY
	 Steps our reduced system forward in time using a reduced model
	 */
	virtual Eigen::VectorXd reduced_step(const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& bc);

	/*
	Steps our full system forward in time,
	u_curr - 3nx1 vector of current displacement
	u_prev - 3nx1 vector of previous displacement
	bc - 3cx1 vector of constrained displacements. Needs to match the number of rows preivously given as input in S

	returns:
	u_next - 3nx1 vector of simulation displacements that satisfy Su_next = bc and minimize our corotational elastic energy + inertia
	*/
	virtual Eigen::VectorXd full_step(const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev, const Eigen::VectorXd& bc);

	/*
	Also recomputes energy and grad_e
	*/
	virtual Eigen::VectorXd reduced_step_eval_metrics(const Eigen::VectorXd& z_curr, const Eigen::VectorXd& z_prev, const Eigen::VectorXd& bc, double& e, double& grad_e);


	/*
	Also returns e and grad_e
	*/
	virtual Eigen::VectorXd full_step_eval_metrics(const Eigen::VectorXd& u_curr, const Eigen::VectorXd& u_prev, const Eigen::VectorXd& bc, double& e, double& grad_e);
	

	virtual void init_modes(int num_modes); //shared with parent

	//slightly more of a priority
	virtual void init_clusters(int num_clusters, int num_feature_modes); //.shared with parent

	virtual void init_system_matrices();
	 
	virtual void init_full_system_matrices();


	virtual void init_reduction_matrices();

	virtual void init_clustering_matrices();

	virtual void full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

	virtual void reduced_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

	virtual void precompute_solvers();


	// These should all be declared here and in none of the children
	//THESE ones are to update the sim in real time... dNot a priority
	virtual void switch_clustering(bool do_clustering);

	virtual void switch_reduction(bool do_reduction);

	virtual void update_material_properties(double ym, double pr);

	virtual void update_equality_constraint(const Eigen::SparseMatrix<double>& S);

	virtual void update_timestep(double dt);

	virtual void update_modes(double new_num_modes);

	virtual void update_clusters(double new_num_clusters);

	virtual void update_modes_cache_dir(std::string new_mode_dir);

	virtual void update_clusters_cache_dir(std::string new_clusters_dir);

	virtual void energy(const Eigen::VectorXd& u, double& bending, double& volume, double& inertia);


public:
	//Initial simulation mesh
	Eigen::MatrixXd X;
	Eigen::MatrixXi T;
	Eigen::MatrixXi F;
	Eigen::VectorXi ext_ind, FiT;

	//Rig jacobian
	Eigen::SparseMatrix<double> J;

	//Rig Null space 3Vx3V that simply sets all rig-controlled entries to 0, and the others are left unchanged.
	Eigen::SparseMatrix<double> N;

	//cache paths
	std::string clusters_file_dir;
	std::string modes_file_dir;

	//reduction parameters
	bool do_reduction;
	int num_modes;
	Eigen::MatrixXd B;

	Eigen::VectorXd L;
	Eigen::MatrixXd B_full;
	Eigen::VectorXd L_full;

	//clustering parameters
	bool do_clustering;
	int num_clusters;
	Eigen::VectorXi labels;


	bool do_inertia;

	//sim parameters
	double stiffness, incompressibility;
	double dt;


	int num_modal_features;

	struct FullSimPrecomputedMatrices
	{
		Eigen::SparseMatrix<double> M; //3|V|x3|V| lumped massmatrix
		Eigen::SparseMatrix<double> A; // 3|V|x3|V|system matrix (stiffness * C + 1/h^2 M)

		Eigen::SparseMatrix<double> Vol; //|T|x|T| diagonal matrix of tetrahedral volumes
		Eigen::SparseMatrix<double> Vol_c; //c x c diagonal matrix of tet volumes
		Eigen::SparseMatrix<double> Vol_exp; //9|T|x9|T| diagonal matrix of expanded tetrahedral volumes 
		
		
		Eigen::SparseMatrix<double> K; //9|T|x3|V|maps displacements contributions to their flattened deformation gradient
		Eigen::VectorXd H;				//constant term required in obtaining defo gradient from u. F = H + K*u;


		Eigen::SparseMatrix<double> KMK; // K^T * Vol_exp * K  -> frequently used quantity, actually should be equivalent to cotan laplacian
		Eigen::VectorXd KMH;			 // K^T * Vol_exp * H

		Eigen::SparseMatrix<double> traceMat; // |T|x|9T| trace matrix operator that takes all flattened 3x3 entries and calculates the trace
												// vec(tr(F)) = traceMat * vec(F)
	
		//Eigen::SparseMatrix<double> G;			// c x |T| grouping matrix that sums along tets in each cluster
		Eigen::SparseMatrix<double> G_1;		// 9c x 9|T|  grouping matrix that sums along  tets in each cluster.. but expanded to work on flattened 3x3 per tet quantities
		Eigen::SparseMatrix<double> G_m;		// grouping matrix that does a mass weighted sum along tets in each cluster. This is expanded to work on flattened 3x3 per cluster quantities
		double HMH;								// H.transpose()*Vol_exp * H
	
		//clusters precomputation
		Eigen::VectorXd GmH;					//G_m * H
		Eigen::SparseMatrix<double> GmK;		// G_m * K

		Eigen::VectorXd GMH;					//G_1 * Vol_exp * H
		Eigen::SparseMatrix<double> GMK;		//G_1 * Vol_exp * K

		Eigen::SparseMatrix<double> KMG;		//K^T * Vol_exp * G_1

		//Following terms are used in the reduced energy calculation
		Eigen::SparseMatrix<double> GMG;		// G_1.transpose() * Vol_exp * G_1
		Eigen::RowVectorXd HMG;					// H.transpose() * Vol_exp * G_1

	}fmp;										


	struct ReducedSimPrecomputedMatrices
	{
		Eigen::MatrixXd BAB;				//B^T A B
		//energy computation
		Eigen::MatrixXd BMB;				// B^T M * B
		Eigen::MatrixXd BKMKB;				// B^T K^T Vol_exp K B
		Eigen::VectorXd BKMH;				// B^T K^T Vol_exp H
		Eigen::MatrixXd BKMG;				// B^T K^T Vol_exp G
		//gradient computation
		Eigen::MatrixXd GmKB;				//G_m^T K B

		Eigen::MatrixXd SB;
	}rmp;

	//quantities that are computed once at the start of each timestep
	struct DynamicMatrices
	{
		Eigen::VectorXd y;
		Eigen::VectorXd BMy;				//B^T M y	
		double yMy;							// y^T M y
	} dm;
	QuasiNewtonSolver* full_newton_solver;
	DenseQuasiNewtonSolver* reduced_newton_solver;

};
