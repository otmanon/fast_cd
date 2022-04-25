#pragma once
#include <Eigen/Core>
#include <iostream>
#include <fstream>

bool write_binary(std::string filename, const Eigen::MatrixXd& A);
bool read_binary(std::string filename, Eigen::MatrixXd& A);