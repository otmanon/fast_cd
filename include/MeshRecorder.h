#pragma once
#include <Eigen/Core>
#include <igl/writeMSH.h>
#include <igl/writeOBJ.h>
#include <filesystem>

class MeshRecorder
{
	public:
		std::vector<Eigen::MatrixXd> V;
		std::vector<Eigen::MatrixXi>  F;

		std::vector<Eigen::MatrixXd> TC;
		std::vector<Eigen::MatrixXi>  FTC;

		std::vector<Eigen::MatrixXd> N;
		std::vector<Eigen::MatrixXi>  FN;
		char name[8] = "";


		MeshRecorder()
		{
		
		}

		void record_frame(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
		{
			this->V.push_back(V);
			this->F.push_back(F);
		}

		void record_frame(MatrixXd& V, MatrixXd& TC, MatrixXd& N, MatrixXi& F, MatrixXi& FTC, MatrixXi& FN)
		{
			this->V.push_back(V);
			this->TC.push_back(TC);
			this->N.push_back(N);
			this->F.push_back(F);
			this->FTC.push_back(FTC);
			this->FN.push_back(FN);
		}

		bool save(std::string dir)
		{
			namespace fs = std::filesystem;
			bool success = true;
			if (!fs::exists(fs::path(dir)))
			{
				fs::create_directories(fs::path(dir));
			}
		
			for (int i = 0; i < this->V.size(); i++)
			{
				sprintf(name, "%04i", i);
				std::string filepath = dir + "/" + name + ".obj";
				Eigen::MatrixXd V2 = V[i];
				Eigen::MatrixXi F2 = F[i];
				success = success && igl::writeOBJ(filepath, V2, F2);
			}
			return success;
		}

		bool save_texture(std::string dir)
		{
			namespace fs = std::filesystem;
			bool success = true;
			if (!fs::exists(fs::path(dir)))
			{
				fs::create_directories(fs::path(dir));
			}

			for (int i = 0; i < this->V.size(); i++)
			{
				sprintf(name, "%04i", i);
				std::string filepath = dir + "/" + name + ".obj";
				Eigen::MatrixXd V2 = V[i];
				Eigen::MatrixXi F2 = F[i];

				Eigen::MatrixXd N2 = N[i];
				Eigen::MatrixXi FN2 = FN[i];

				Eigen::MatrixXd TC2 = TC[i];
				Eigen::MatrixXi FTC2 = FTC[i];
				success = success && igl::writeOBJ(filepath, V2, F2, N2, FN2, TC2, FTC2);

			}
			return success;
		}
};