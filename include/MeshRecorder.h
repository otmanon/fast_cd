#pragma once
#include <Eigen/Core>
#include <igl/writeMSH.h>
#include <igl/writeOBJ.h>
#include <filesystem>

class MeshRecorder
{
public:
	std::vector<Eigen::MatrixXd> V;
	std::vector<Eigen::MatrixXi> T, F;

	MeshRecorder()
	{
		
	}

	void record_frame(Eigen::MatrixXd& V, Eigen::MatrixXi& T, Eigen::MatrixXi& F)
	{
		this->V.push_back(V);
		this->T.push_back(T);
		this->F.push_back(F);
	}

	void save(std::string dir)
	{
		namespace fs = std::filesystem;

		if (!fs::exists(fs::path(dir)))
		{
			fs::create_directories(fs::path(dir));
		}
		char name[8] = "";
		for (int i = 0; i < this->V.size()-1; i++)
		{
			sprintf(name, "%04i.png", i);
			std::string filepath = dir + "/" + name + ".obj";
			Eigen::MatrixXd V2 = V[i];
			Eigen::MatrixXi F2 = F[i];
			igl::writeOBJ(filepath, V2, F2);
		}
	}
};