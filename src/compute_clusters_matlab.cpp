
#include "compute_clusters_matlab.h"
#include <iostream>
#include <cassert>
#include <igl/sort.h>
#include <igl/get_seconds.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/average_onto_faces.h>


void compute_clusters_matlab(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C)
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

   // igl::kmeans(B_faces, num_clusters, C, labels);
  
    
    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine, "B", B_faces);
   // igl::matlab::mlsetmatrix(&engine, "B", B);
    igl::matlab::mlsetscalar(&engine, "k", num_clusters);
    igl::matlab::mleval(&engine, "tic");
    igl::matlab::mleval(&engine, "[idx, C] = kmeans(B,k)");
    igl::matlab::mleval(&engine, "I = real(idx)");
    igl::matlab::mleval(&engine, "C = real(C)");
    igl::matlab::mleval(&engine, "time = toc");
    // igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
    Eigen::MatrixXi labels_mat;
    igl::matlab::mlgetmatrix(&engine, "I", labels_mat);
    igl::matlab::mlgetmatrix(&engine, "C", C);
    labels = labels_mat.col(0);//
    Eigen::MatrixXd S_mat;
    double t;
    igl::matlab::mlgetmatrix(&engine, "S", S_mat);
}


void compute_clusters_matlab_initguess(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::MatrixXd& L, const int num_clusters, const int num_feature_modes, const Eigen::MatrixXd& C0, Eigen::VectorXi& labels, Eigen::MatrixXd& C)
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

    // igl::kmeans(B_faces, num_clusters, C, labels);


    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine, "B", B_faces);
    igl::matlab::mlsetmatrix(&engine, "C0", C0);
    // igl::matlab::mlsetmatrix(&engine, "B", B);
    igl::matlab::mlsetscalar(&engine, "k", num_clusters);
    igl::matlab::mleval(&engine, "tic");
    igl::matlab::mleval(&engine, "[idx, C] = kmeans(B,k, 'Start', C0)");
    igl::matlab::mleval(&engine, "I = real(idx)");
    igl::matlab::mleval(&engine, "C = real(C)");
    igl::matlab::mleval(&engine, "time = toc");
    // igl::matlab::mleval(&engine, "save(\"aca_modes.mat\")");
    Eigen::MatrixXi labels_mat;
    igl::matlab::mlgetmatrix(&engine, "I", labels_mat);
    igl::matlab::mlgetmatrix(&engine, "C", C);
    labels = labels_mat.col(0);//
    Eigen::MatrixXd S_mat;
    double t;
    igl::matlab::mlgetmatrix(&engine, "S", S_mat);
}

