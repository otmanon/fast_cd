#pragma once
#include <Eigen/Core>
#include <igl/writeDMAT.h>
#include <filesystem>
class ScalarRecorder
{//
public:
	Eigen::VectorXd P;


	ScalarRecorder()
	{
		
	}

	void record_frame(double n)
	{
		P.conservativeResize(P.rows()+ 1);
		P(P.rows() - 1) = n;
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

	double mean()
	{
		if (P.rows() > 0)
			return P.mean();
		else
			return 0.0;
	}

	/*
	Returns mean only over last n frames
	*/
	double mean(int n)
	{
		if (P.rows() > n)
		{
			VectorXd p = P.bottomRows(n);
			return p.mean();
		}
		else
			return 0.0;
	}
};