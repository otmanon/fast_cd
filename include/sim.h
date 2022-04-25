#pragma once
#include <Eigen\Core>
#include <Eigen\Sparse>

#include "CDNewtonSolver.h"
#include "FastCDNewtonSolver.h"
//TODO: Make THIS the parent class of fast simulators, and have FastCD Sim inherit from it
class Sim
{
public:
	//TODO: need to  implement all these method
	Sim();

	/*
	Creats
	X -  nx3 rest mesh pose
	T -  |T|x4 tet indices into X
	S -  (c x 3n) linear equality constraints on flattened X such that Su = f specifies vertex displacements of f for all constrained vertices

	ym - youngs modulus
	pr - poisson ratio
	dt - timestep
	*/
	Sim(const Eigen::MatrixXd& X,const  Eigen::MatrixXi& T,const  Eigen::SparseMatrix<double>& S, double ym, double pr, double dt);

	void step(const Eigen::VectorXd& u_curr,const  Eigen::VectorXd& u_prev, const Eigen::VectorXd& bc);

	void init_system_matrices(); //overrides parent

	void init_full_system_matrices();

	void full_qnewton_energy_grad(std::function<double(const Eigen::VectorXd&)>& f, std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

	void precompute_solvers();

	//THESE ones are to update the sim in real time... dNot a priority
	void switch_clustering(bool do_clustering);

	void switch_reduction(bool do_reduction);

	void update_material_properties(double ym, double pr);

	void update_timestep(double dt);

	void update_equality_constraints(const Eigen::SparseMatrix<double>& S, std::string new_modes_dir, std::string new_clusters_dir);

	void update_modes_cache_dir(std::string new_mode_dir);

	void update_clusters_cache_dir(std::string new_clusters_dir);

public:
	//Initial simulation mesh
	Eigen::MatrixXd X;
	Eigen::MatrixXi T;
	Eigen::MatrixXi F;

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

	Eigen::VectorXd S;
	Eigen::MatrixXd B_full;
	Eigen::VectorXd L_full;


	//sim parameters
	double stiffness, incompressibility;
	double dt;

	struct FullSimPrecomputedMatrices
	{

	}fmp;



	CDNewtonSolver* newton_solver;

};
