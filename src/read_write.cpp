#include "read_write.h"
using namespace Eigen;
bool write_binary(std::string filename, const MatrixXd& A) {
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    typename MatrixXd::Index rows = A.rows(), cols = A.cols();
    out.write((char*)(&rows), sizeof(typename MatrixXd::Index));
    out.write((char*)(&cols), sizeof(typename MatrixXd::Index));
    out.write((char*)A.data(), rows * cols * sizeof(typename MatrixXd::Scalar));
    out.close();
    return true;
}


bool read_binary(std::string filename, MatrixXd& A) {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in.fail()) return false;
    typename MatrixXd::Index rows = 0, cols = 0;
    in.read((char*)(&rows), sizeof(typename MatrixXd::Index));
    in.read((char*)(&cols), sizeof(typename MatrixXd::Index));
    A.resize(rows, cols);
    in.read((char*)A.data(), rows * cols * sizeof(typename MatrixXd::Scalar));
    in.close();
    return true;
}