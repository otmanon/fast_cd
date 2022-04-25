#pragma once
#include <Eigen/core>
#include <iostream>
#include <fstream>

bool write_binary(std::string filename, const Eigen::MatrixXd& A);
bool read_binary(std::string filename, Eigen::MatrixXd& A);