#pragma once
#include <Eigen/Core>
#include <igl/writeDMAT.h>
#include <filesystem>
class VectorRecorder
{//
public:
	Eigen::MatrixXd P;

	VectorRecorder() {};
	VectorRecorder(Eigen::VectorXd& p0)
	{
		P.resize(p0.rows(), 0);
	}

	VectorRecorder(int n)
	{
		P.resize(n, 0);
	}

	void record_frame(Eigen::VectorXd& p)
	{
		P.conservativeResize(P.rows(), P.cols() + 1);
		P.col(P.cols() - 1) = p;
	}

	void save(std::string filepath)
	{
		namespace fs = std::filesystem;

		if (!fs::exists(fs::path(filepath).parent_path()))
		{
			fs::create_directories(fs::path(filepath).parent_path());
		}
		igl::writeDMAT(filepath, P);
	}
};