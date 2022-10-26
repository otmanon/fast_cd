#pragma once
#include <Eigen/Core>
#include <string>
#include <igl/readDMAT.h>

using namespace std; using namespace Eigen;

/*
From cache_dir, reads all matrices and vectors required for the reduced static precomputation of sp
*/
bool read_fast_cd_sim_static_precomputation(string& cache_dir, MatrixXd& B, VectorXd& L,  VectorXi& l, 
	MatrixXd& BCB, MatrixXd& BMB, MatrixXd& BAB, MatrixXd& AeqB, MatrixXd& GmKB, MatrixXd& GmKJ, 
	VectorXd& GmKx, MatrixXd& G1VKB, MatrixXd& BMJ, VectorXd& BMx,MatrixXd& BCJ, VectorXd& BCx )
{
	bool t = true;
	//cache
	t = igl::readDMAT(cache_dir + "/B.DMAT", B) && t;
	t = igl::readDMAT(cache_dir + "/L.DMAT", L) && t;
	VectorXd labels;
	t = igl::readDMAT(cache_dir + "/labels.DMAT", labels) && t;
	l = labels.cast<int>();
	t = igl::readDMAT(cache_dir + "/BCB.DMAT", BCB) && t;
	t = igl::readDMAT(cache_dir + "/BMB.DMAT", BMB) && t;
	t = igl::readDMAT(cache_dir + "/BAB.DMAT", BAB) && t;
	t = igl::readDMAT(cache_dir + "/AeqB.DMAT", AeqB) && t;
	t = igl::readDMAT(cache_dir + "/GmKB.DMAT", GmKB) && t;
	t = igl::readDMAT(cache_dir + "/GmKJ.DMAT", GmKJ) && t;
	t = igl::readDMAT(cache_dir + "/GmKx.DMAT", GmKx) && t;
	t = igl::readDMAT(cache_dir + "/G1VKB.DMAT", G1VKB) && t;
	t = igl::readDMAT(cache_dir + "/BMJ.DMAT", BMJ) && t;
	t = igl::readDMAT(cache_dir + "/BMx.DMAT", BMx) && t;
	t = igl::readDMAT(cache_dir + "/BCJ.DMAT", BCJ) && t;
	t = igl::readDMAT(cache_dir + "/BCx.DMAT", BCx) && t;
	return t;
}