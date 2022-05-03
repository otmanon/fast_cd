#include "compute_clusters_igl.h"
#include "kmeans.h"
#include <igl/average_onto_faces.h>

void compute_clusters_igl(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C)
{
    int c, k;
    Eigen::MatrixXd B_tmps;
    c = num_feature_modes < B.cols() ? num_feature_modes : B.cols();
    Eigen::RowVectorXd L_tmp = L.topRows(c).transpose();
    Eigen::MatrixXd B_tmp = B.block(0, 0, B.rows(), num_feature_modes);
    B_tmp.array().rowwise() /= L_tmp.row(0).array();
    B_tmp.array().rowwise() /= L_tmp.row(0).array();
    Eigen::MatrixXd B_verts = Eigen::Map<Eigen::MatrixXd>(B_tmp.data(), B_tmp.rows() / 3, B_tmp.cols() * 3);
    //optionally normalize by S
    Eigen::MatrixXd B_faces;
    igl::average_onto_faces(T, B_verts, B_faces);
    igl::kmeans(B_faces, num_clusters, C, labels);
}