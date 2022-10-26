#pragma once
#include <Eigen/Core>
#include <igl/writeDMAT.h>
#include <filesystem>
using namespace std; using namespace Eigen;

/*
From cache_dir, reads all matrices and vectors required for the reduced static precomputation of sp
*/
bool write_fast_cd_sim_static_precomputation(string& cache_dir, MatrixXd& B, VectorXd& L, VectorXi& l,
	MatrixXd& BCB, MatrixXd& BMB, MatrixXd& BAB, MatrixXd& AeqB, MatrixXd& GmKB, MatrixXd& GmKJ,
	VectorXd& GmKx, MatrixXd& G1VKB, MatrixXd& BMJ, VectorXd& BMx, MatrixXd& BCJ, VectorXd& BCx)
{
	namespace fs = std::filesystem;

	if (!fs::exists(fs::path(cache_dir)))
		fs::create_directories(fs::path(cache_dir));
	bool t = true;
	//cache
	t = igl::writeDMAT(cache_dir + "/B.DMAT", B, false) && t;
	t = igl::writeDMAT(cache_dir + "/L.DMAT", L,  false) && t;
	VectorXd labels = l.cast<double>();
	t = igl::writeDMAT(cache_dir + "/labels.DMAT", labels, false) && t;

	t = igl::writeDMAT(cache_dir + "/BCB.DMAT", BCB, false) && t;
	t = igl::writeDMAT(cache_dir + "/BMB.DMAT", BMB, false) && t;
	t = igl::writeDMAT(cache_dir + "/BAB.DMAT", BAB, false) && t;
	t = igl::writeDMAT(cache_dir + "/AeqB.DMAT", AeqB, false) && t;
	t = igl::writeDMAT(cache_dir + "/GmKB.DMAT", GmKB, false) && t;
	t = igl::writeDMAT(cache_dir + "/GmKJ.DMAT", GmKJ, false) && t;
	t = igl::writeDMAT(cache_dir + "/GmKx.DMAT", GmKx, false) && t;
	t = igl::writeDMAT(cache_dir + "/G1VKB.DMAT", G1VKB, false) && t;
	t = igl::writeDMAT(cache_dir + "/BMJ.DMAT", BMJ, false) && t;
	t = igl::writeDMAT(cache_dir + "/BMx.DMAT", BMx, false) && t;
	t = igl::writeDMAT(cache_dir + "/BCJ.DMAT", BCJ, false) && t;
	t = igl::writeDMAT(cache_dir + "/BCx.DMAT", BCx, false) && t;
	return t;
}