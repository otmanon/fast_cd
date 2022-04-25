#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Eigen;
/*
Given system Ax=b, enforce linear equality constraints C^Tx = d through lagrange multipliers. This involves constructing
the KKT system:
[ A   CT] [x] = [b]
[ C    0] [l] = [d]

Which becomes
Q y = f

*/
void augment_with_linear_constraints(const SparseMatrix<double>& A, VectorXd& b, 
	const SparseMatrix<double>& Aeq, const VectorXd& beq, const SparseMatrix<double>& Q, VectorXd& d);

/*
Given system Ax=b, enforce linear equality constraints C^Tx = d through lagrange multipliers. This involves constructing
the KKT system:
[ A   CT] [x] = [b]
[ C    0] [l] = [d]

Which becomes
Q y = f

This method only returns the reduced system matrix. 
*/
void augment_with_linear_constraints(const SparseMatrix<double>& A, const SparseMatrix<double>& Aeq, SparseMatrix<double>& Q);


void augment_with_linear_constraints(const Eigen::MatrixXd& A, const  Eigen::MatrixXd& C, Eigen::MatrixXd& Q);