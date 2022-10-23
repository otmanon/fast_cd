#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

struct arap_fast_cd_sim
{

};

struct arap_fast_cd_sim_state
{

};

struct arap_fast_cd_sim_params
{
	//timestep
	double h;
	//youngs modulus (one per each tet)
	Eigen::VectorXd ym;
	//poisson's ratio
	Eigen::VectorXd pr;

};

struct local_global_solver_params
{

};

struct arap_fast_cd_solver
{

};


struct arap_fast_cd_static_precomp
{

};


struct arap_fast_cd_dynamic_precomp
{

};