#include "HandleRig.h"
#include "interweaving_matrix.h"
#include <igl/slice.h>
#include <igl/boundary_facets.h>
#include <igl/centroid.h>
#include <igl/unique.h>

#include <filesystem>

#include "compute_handle_positions_from_rig_parameters.h"
#include <igl/repdiag.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
#include <igl/list_to_matrix.h>
#include <igl/matrix_to_list.h>
HandleRig::HandleRig(std::string filename, double radius)
{
	read_rig_from_json(filename);
	init_rig_jacobian();
	init_null_space();

	init_rig_selection_matrix(radius);
}

HandleRig::HandleRig(Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius)
{
	this->X = X;
	this->W = SparseMatrix<double>(X.rows(), 1);
	this->W.setConstant(1);

	Eigen::MatrixXi F; Eigen::VectorXi n1, n2;
	igl::boundary_facets(T, F, n1, n2);

	Eigen::Vector3d centroid;
	igl::centroid(X, F, centroid);
	Eigen::Matrix4f M;
	M.setIdentity();
	M.block(0, 3, 3, 1) = centroid.cast<float>();

	//need to row order flatten this.
	Eigen::MatrixXf A = M.topRows(3);
	A = A.transpose().eval();

	p0 = Eigen::Map<Eigen::VectorXf>(A.data(), 12).cast<double>();
	init_rig_jacobian();
	init_null_space();

	rig_pinning = "ball";
	init_rig_selection_matrix(radius);


}

HandleRig::HandleRig(Eigen::MatrixXd& X, Eigen::VectorXd& p, Eigen::MatrixXd& W, double radius)
{
	this->p0 = p;
	this->X = X;
	this->W = W;
	init_rig_jacobian();
	init_null_space();

	init_rig_selection_matrix(radius);
	//get exterior indices
	//Eigen::VectorXi bI;
	//Eigen::MatrixXi F;
	////igl::boundary_facets(T, F);
	//igl::unique(F, bI);
	//also keep track of exterior slices of jacobian
}

/*
Finds indices of vertex that are bound by the rig. A vertex is bound by the rig if it is within a radius of r from the origin, where r is taken as mesh_width*0.05
*/
void HandleRig::init_rig_selection_matrix(double radius)
{

	bI.resize(W.cols()); // get one index list per bone

	int ci = 0; //keeps track of the number of constrained vertices
	if (rig_pinning == "ball")
	{
		double r = (X.colwise().maxCoeff() - X.colwise().minCoeff()).maxCoeff() * radius;
		Eigen::MatrixXd V;
		compute_handle_positions_from_parameters(p0, V);
		Eigen::VectorXi VI;
		igl::colon(0, V.rows() - 1, VI);
		//get closest point
		Eigen::VectorXd sqrD;
		Eigen::MatrixXd CP;
		Eigen::VectorXi I;
		igl::point_mesh_squared_distance(X, V, VI, sqrD, I, CP);
		for (int i = 0; i < I.rows(); i++)
		{
		//if closer than r away
			if (sqrD(i) < r * r)
			{
				bI[I(i)].conservativeResize(bI[I(i)].size() + 1);
				bI[I(i)](bI[I(i)].size() - 1) = i;
				ci += 1;
			}
		}
	}if (rig_pinning == "weight_threshold")
	{
		int ci = 0; //keeps track of the number of constrained vertices
		for (int b = 0; b < W.cols(); b++)
		{
			for (int i = 0; i < W.rows(); i++)
			{
				//if closer than r away
				if (W(i, b) > 0.8)  //there will not be multiple rows that have this happen
				{
					bI[b].conservativeResize(bI[b].size() + 1);
					bI[b](bI[b].size() - 1) = i;
					ci += 1;
				}
			}
		}
	}
	if (rig_pinning == "slab")
	{
		Eigen::MatrixXd V;
		compute_handle_positions_from_parameters(p0, V);
		Eigen::VectorXi VI;
		igl::colon(0, V.rows() - 1, VI);
		//get closest point
		Eigen::VectorXd sqrD;
		Eigen::MatrixXd CP;
		Eigen::VectorXi I;
		igl::point_mesh_squared_distance(X, V, VI, sqrD, I, CP);
		for (int i = 0; i < I.rows(); i++)
		{
			//if closer than r away
			if (abs(X(i, 0) - V(I(i), 0)) < 1e-1)
			{
				bI[I(i)].conservativeResize(bI[I(i)].size() + 1);
				bI[I(i)](bI[I(i)].size() - 1) = i;
				ci += 1;
			}
		}
	}
	/*Can also do this by weight, or by slab One approach, there are a few others we can try

	}
	*/
	std::vector<Eigen::Triplet<double>> tripletList;
	int count = 0;
	//go through list of bones
	for (int b = 0; b < bI.size(); b++)
	{
		//go through list of indices constrained by bones
		for (int i = 0; i < bI[b].rows(); i++)
		{
			tripletList.push_back(Eigen::Triplet<double>(i + 0 * ci + count, bI[b](i) + 0 * X.rows(), 1));
			tripletList.push_back(Eigen::Triplet<double>(i + 1 * ci + count, bI[b](i) + 1 * X.rows(), 1));
			tripletList.push_back(Eigen::Triplet<double>(i + 2 * ci + count, bI[b](i) + 2 * X.rows(), 1));
	
		}
		count += bI[b].rows();
	}

	S.resize(3 * ci, 3 * X.rows());
	S.setFromTriplets(tripletList.begin(), tripletList.end());

	Eigen::SparseMatrix<double> r2c;
	interweaving_matrix(S.rows()/3, 3, r2c);

	//S = (r2c.transpose()) * S;
}


HandleRig::HandleRig(Eigen::MatrixXd& X,  std::vector<Eigen::Matrix4f>& P, Eigen::MatrixXd& W, double radius)
{
	this->X = X;
	this->W = W;
	Eigen::VectorXd p;
	p.resize(P.size() * 12);
	int rows = 3;
	int cols = 4;
	for (int i = 0; i < P.size(); i++)
	{
		//do row-order flattening of rig parameters (as Jacobian expects)
		const Eigen::Matrix4d Q = P[i].cast<double>();
		for (int row = 0; row < Q.rows()-1; row++)
		{
			p.block(i * cols * rows + row * cols, 0, cols, 1) = Q.row(row).transpose();
		}
	}

	p0 = p;
	init_rig_jacobian();
	init_null_space();
	init_rig_selection_matrix(radius);
	//get exterior indices

	//also keep track of exterior slices of jacobian
	//igl::slice(J, bI, 1, J_ext);
}

void HandleRig::init_rig_jacobian()
{
	double w_ij;
	Eigen::Vector3d vi, wvi;
	int b;
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(X.rows() *12* W.cols()); // 12 entries for each bone-vertex pair
	int v_dim = 3;
	int a_dim = 4;
	int row, col;

	for (int i = 0; i < X.rows(); i++)
	{
		vi = X.row(i);
		for (int b = 0; b < W.cols(); b++)
		{
			w_ij = W(i, b);
			wvi = w_ij * vi;
			//each bone-vertex pair contributes to a 3x12 block in the jacobian
			for (int ci = 0; ci < v_dim; ci++)
			{
				row = i * v_dim + ci;
				for (int cj = 0; cj < a_dim - 1; cj++)
				{
					col = b * a_dim * v_dim + cj + ci * a_dim;
					tripletList.emplace_back(Eigen::Triplet<double>(row, col, wvi(cj)));
				}
				col = b * a_dim * v_dim + a_dim - 1 + ci * a_dim;
				tripletList.emplace_back(Eigen::Triplet<double>(row, col, w_ij));
			}
		}
	}
	J.resize(v_dim * X.rows(), W.cols()*12);
	J.setFromTriplets(tripletList.begin(), tripletList.end());

	Eigen::SparseMatrix<double> S;
	interweaving_matrix(X.rows(), X.cols(), S);
	J = S * J;
}


void HandleRig::get_rig_jacobian(Eigen::SparseMatrix<double>& J)
{
	J = this->J;
}


void HandleRig::get_rig_motion(Eigen::VectorXd& p, Eigen::VectorXd& ur)
{
	ur = J * p ;
}



bool HandleRig::write_rig_to_json(std::string filename)
{
	
	namespace fs = std::filesystem;
	
	filename = fs::path(filename).extension().string() == ".json" ? filename : filename + ".json";

	//fs::path p = fs::path(filename);
	
	//std::string parent_path = p.parent_path().string();
	if (!fs::exists(fs::path(filename).parent_path()))
	{
		
		fs::create_directories(fs::path(filename).parent_path());
	}

	//write rig to json
	std::ofstream o(filename);
	json j;

	std::vector<std::vector<double>> X_list, W_list;
	std::vector<double> p0_list;
        igl::matrix_to_list(X,X_list);
	igl::matrix_to_list(W,W_list);
	p0_list = std::vector<double>(p0.data(), p0.data() + p0.rows());

	j["rig_type"] = "handle_rig";
	//save mesh positions
	j["X"] = X_list;
	j["W"] = W_list;
	j["p0"] = p0_list;

	o << std::setw(4) << j << std::endl;

	return true;
}


bool HandleRig::read_rig_from_json(std::string filename)
{
	namespace fs = std::filesystem;


	json j;

	std::ifstream i(filename);
	i >> j;

	rig_pinning = j.value("rig_pinning", "ball");
	std::vector<std::vector<double>> X_list = j["X"];
	std::vector<std::vector<double>> W_list = j["W"];
	std::vector<double> p0_list = j["p0"];

	igl::list_to_matrix(X_list,X);
	igl::list_to_matrix(W_list,W);
	p0 = Eigen::Map<Eigen::VectorXd>(&p0_list[0], p0_list.size());
	return true;
}

void HandleRig::init_null_space()
{
	Eigen::SparseMatrix<double> N_small;
	N_small.resize(X.rows(), X.rows());
	std::vector<Eigen::Triplet<double>> tripletList;
	//tripletList.reserve();
	for (int i = 0; i < W.rows(); i++)
	{
		double w = W.row(i).sum();
		if (w < 1e-8)
		{
			//append to diagonal
			tripletList.push_back(Eigen::Triplet<double>(i, i, 1.0));
		}
	}
	N_small.setFromTriplets(tripletList.begin(), tripletList.end());

	N = igl::repdiag(N_small, 3);

}
