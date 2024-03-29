#include "compute_clusters_igl.h"
#include "kmeans.h"
#include <igl/average_onto_faces.h>
#include "split_components.h"
#include "tet_adjacency_matrix.h"
void compute_clusters_igl(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C, bool split_components)
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
    Eigen::VectorXi labels_global;
    igl::kmeans(B_faces, num_clusters, C, labels_global);

    Eigen::SparseMatrix<bool, 0, int> A;
    tet_adjacency_matrix_vert_share(T, A); 
    if (split_components)
        igl::split_components(labels_global, A, labels);
        else
        labels = labels_global;
}



void compute_clusters_displacement_features(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C, bool split_components)
{
    int c, k;
    Eigen::MatrixXd B_tmps;
    c = num_feature_modes < B.cols() ? num_feature_modes : B.cols();
    Eigen::RowVectorXd L_tmp = L.topRows(c).transpose();
    Eigen::MatrixXd B_tmp = B.block(0, 0, B.rows(), num_feature_modes);
    B_tmp.array().rowwise() /= L_tmp.row(0).array().abs();
    B_tmp.array().rowwise() /= L_tmp.row(0).array().abs();
    Eigen::MatrixXd B_verts = Eigen::Map<Eigen::MatrixXd>(B_tmp.data(), B_tmp.rows() / 3, B_tmp.cols() * 3);
    //optionally normalize by S
    Eigen::MatrixXd B_faces;
    igl::average_onto_faces(T, B_verts, B_faces);
    Eigen::VectorXi labels_global;
    igl::kmeans(B_faces, num_clusters, C, labels_global);

    Eigen::SparseMatrix<bool, 0, int> A;
    tet_adjacency_matrix_vert_share(T, A);

    if (split_components)
        igl::split_components(labels_global, A, labels);
    else
        labels = labels_global;
}


void compute_clusters_weight_features(const Eigen::MatrixXi& T, const Eigen::MatrixXd& W, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C, bool split_components)
{
    int c, k;
    Eigen::MatrixXd B_tmps;
    c = num_feature_modes < W.cols() ? num_feature_modes : W.cols();
    Eigen::RowVectorXd L_tmp = L.topRows(c).transpose();
    Eigen::MatrixXd W_tmp = W.block(0, 0, W.rows(), num_feature_modes);
    W_tmp.array().rowwise() /= L_tmp.row(0).array();
    W_tmp.array().rowwise() /= L_tmp.row(0).array();
    //optionally normalize by S
    Eigen::MatrixXd B_faces;
    igl::average_onto_faces(T, W_tmp, B_faces);
    Eigen::VectorXi labels_global;
    igl::kmeans(B_faces, num_clusters, C, labels_global);

    Eigen::SparseMatrix<bool, 0, int> A;
    tet_adjacency_matrix_vert_share(T, A);
    
    if (split_components)
        igl::split_components(labels_global, A, labels);
    else
        labels = labels_global;
}