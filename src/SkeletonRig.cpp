#include "SkeletonRig.h"
#include "update_parameters_at_handle.h"
#include "surface_to_volume_weights.h"
#include "matrix4f_from_parameters.h"

#include <igl/boundary_facets.h>
#include <igl/procrustes.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include "prolongation.h"
#include <igl/find.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice_into.h>
#include <igl/list_to_matrix.h>
#include <igl/colon.h>
#include <igl/matrix_to_list.h>
#include "compute_handle_positions_from_rig_parameters.h"
#include "get_tip_positions_from_parameters.h"
#include "get_joint_positions_from_parameters.h"
#include "interweaving_matrix.h"
#include <filesystem>
#include "igl/writeDMAT.h"
#include <Eigen/Geometry>
#include <igl/writeMSH.h>

void SkeletonRig::init_rig_selection_matrix(double radius)
{

	bI.resize(W.cols()); // get one index list per bone

	int ci = 0; //keeps track of the number of constrained vertices

		Eigen::MatrixXd joints, tips, C;
		get_joint_positions_from_parameters(p0, joints);
		get_tip_positions_from_parameters(p0, lengths, tips);
		C.resize(joints.rows() + tips.rows(), 3);
		C.topRows(joints.rows()) = joints;
		C.bottomRows(tips.rows()) = tips;

		Eigen::MatrixXi BE;
		Eigen::VectorXi BE1, BE2;
		igl::colon<int>(0, joints.rows() - 1, BE1);
		igl::colon<int>(joints.rows(), tips.rows() + joints.rows() -1, BE2);
		BE.resize(BE1.rows(), 2);
		BE.col(0) = BE1;
		BE.col(1) = BE2;


		//get closest point
		Eigen::VectorXd sqrD;
		Eigen::MatrixXd CP;
		Eigen::VectorXi I;
		igl::point_mesh_squared_distance(X, C, BE, sqrD, I, CP);
		for (int i = 0; i < I.rows(); i++)
		{
			//if closer than r away
			if (sqrD(i) < radius*radius)
			{
				bI[I(i)].conservativeResize(bI[I(i)].size() + 1);
				bI[I(i)](bI[I(i)].size() - 1) = i;
				ci += 1;
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
	interweaving_matrix(S.rows() / 3, 3, r2c);

	//S = (r2c.transpose()) * S;
}

void read_bones_from_json(json& j, int num_X, Eigen::MatrixXd& surfaceW, Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::VectorXd& lengths)
{
	const int num_bones = j["bones"].size();
	p0.resize(12 * num_bones);
	pI.resize(num_bones); //vector containing parent Indices
	lengths.resize(num_bones); //vector containing edge lengths.
	//use this to convert from flattened matrix to 4x4 affine matrix
	auto parse_affine = [](const json& j) -> Eigen::MatrixX4f
	{
		Eigen::Matrix4f mat;
		mat <<
			j[0][0], j[0][1], j[0][2], j[0][3],
			j[1][0], j[1][1], j[1][2], j[1][3],
			j[2][0], j[2][1], j[2][2], j[2][3],
			0, 0, 0, 1;
		return mat;
	};

	// parse bones to initialize skeleton
	Eigen::VectorXi wI;
	Eigen::VectorXd wV;
	surfaceW.resize(num_X, num_bones);
	surfaceW.setZero();
	int bone_index = 0;
	for (const json& jbone : j["bones"])
	{
		//bone_index = jbone["bone_idx"];  //maybe the bones are ordered in a different way in parent idx than what we loop throughh... so we need this
		pI(bone_index) = jbone["parent_idx"];
		lengths(bone_index) = jbone["bone_length"];
		Eigen::Matrix4f A = parse_affine(jbone["bone_transform"]).matrix();

		update_parameters_at_handle(p0, A, bone_index);
		
		std::vector<int> w_i = jbone["vert_indeces_of_bone"];
		wI = Eigen::Map<Eigen::VectorXi>(w_i.data(), w_i.size());
		std::vector<double> w_v = jbone["weights_of_bone"];
		wV = Eigen::Map<Eigen::VectorXd>(w_v.data(), w_v.size());
		
		Eigen::MatrixXd wV_mat;
		wV_mat = Eigen::Map<Eigen::MatrixXd>(wV.data(), wV.rows(), 1);
		Eigen::VectorXi col(1);
		col(0) = bone_index;
		igl::slice_into(wV_mat, wI, col, surfaceW);
	
		bone_index += 1;
		if (wI.rows() == 0 || wV.rows() == 0)
		{
			std::cout << "Bone " + std::to_string(bone_index) + " has no effect on any vertices. This will result in singular rig Jacobian. Change your rig" << std::endl;
			exit(0);
		}
	}
}

SkeletonRig::SkeletonRig(Eigen::MatrixXd& X, Eigen::VectorXd& p0, Eigen::MatrixXd& W, Eigen::VectorXi& pI, Eigen::VectorXd& lengths, double radius)
{
	this->X = X;
	this->p0 = p0;
	this->W = W;
	this->pI = pI;
	this->lengths = lengths;
	this->radius = radius;
	this->rig_type = "skeleton";
	init_rig_jacobian();
	init_rig_selection_matrix(radius);
}

SkeletonRig::SkeletonRig(std::string surface_file_name, Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius)
{

	printf("Converting a surface rig file format to a volume one. Will attempt to map from surface volume mesh X provided, to the coarse tet mesh we are trying to rig."
				"Then, we will diffuse surface weight values to volume weight values via a diffusion...\n");

		//open surface_file as json
	namespace fs = std::filesystem;

	json j;
	std::ifstream i(surface_file_name);
	i >> j;


	rig_pinning = j.value("rig_pinning", "ball");  //if htis is ball then we treat them like cylinders
	std::vector<std::vector<double>> X_list = j["vertices"];
	Eigen::MatrixXd surface_X;
        igl::list_to_matrix(X_list,surface_X); // hopefully this works
	Eigen::MatrixXd surface_W;
	read_bones_from_json(j, surface_X.rows(), surface_W, this->p0, this->pI, this->lengths);
	this->rig_type = "skeleton";
	//need to do a registration problem between the two methods, for robustness sake... It could be that one is a rotated version of the other.
	Eigen::MatrixXd Xb;
	Eigen::MatrixXi F;
	Eigen::VectorXi n_, FiT, fi;
	igl::boundary_facets(T, F, FiT, n_);
	igl::unique(F, fi);
	igl::slice(X, fi, 1, Xb);

	double scale;      //scale
	Eigen::MatrixXd R; //rotaiton
	Eigen::VectorXd t; //translation
	igl::procrustes(surface_X, Xb, true, false, scale, R, t);
	R *= scale;
	Eigen::MatrixXd Xprime = (surface_X * R).rowwise() + t.transpose();
	int num_b = this->p0.rows() / 12;

	Eigen::Matrix4f A0, A2; Eigen::Matrix3f R0, G;
	Eigen::Affine3f id; id.setIdentity();
	Eigen::Matrix4f Rx =(Eigen::AngleAxisf(-90.0 * igl::PI / 180.0, Eigen::Vector3f::UnitX())*id).matrix();
	Eigen::Matrix4f Rz =  (Eigen::AngleAxisf(90.0 * igl::PI / 180.0, Eigen::Vector3f::UnitZ()) * id).matrix();  //TODO: should have a flat for which axis is the bone aligned on. In the case of blender its on the y axis, but for us its on the x
	// Now we have a rotation, translation and scale... just gotta fit that for all our rig parameters
	for (int i = 0; i < num_b; i++)
	{
		A0 = matrix4f_from_parameters(this->p0, i); //4x4 transformation matrix... now should we rotate first and then translate /scale or what?
		A0 = scale * Rx * A0 * Rz;
		A0.block(0, 3, 3, 1) = A0.block(0, 3, 3, 1) + t.cast<float>();
		update_parameters_at_handle(this->p0, A0, i);
		 A2 = matrix4f_from_parameters(this->p0, i);
	}
	Eigen::SparseMatrix<double> P, PI;
	prolongation(Xprime, X, T, P, PI);

	//finally I can map input W to our "coarse" W that we expect.
	Eigen::MatrixXd W_tmp = P.transpose() * surface_W; //W_tmp will now have values of W on rows mapping to surface_W, and 0 on other rows.
	Eigen::VectorXd sums = P.transpose() * Eigen::VectorXd::Ones(surface_W.rows());
	Eigen::VectorXi all, aI, bI, bJ, bV;
	for (int i = 0; i < sums.rows(); i++)
	{
		if (sums(i) > 1-1e-1)
		{
			sums(i) = 0;
			bI.conservativeResize(bI.rows() + 1);
			bI(bI.size() - 1) = i; //only get indices of non zero entries
		}
	}
	//Eigen::VectorXd W_tmp
	//igl::find(W_tmp, bI, bJ, bV);					//find all non-zero entries here, as they will now be treated as BOUNDARY conditions
	
	//igl::writeDMAT("../data/scene_data/glob_tree/rigs/skeleton_rig/bI.DMAT",  bI);
	//igl::writeDMAT("../data/scene_data/glob_tree/rigs/skeleton_rig/W_tmp.DMAT", W_tmp);
	igl::slice(W_tmp, bI, 1, surface_W);

	//igl::writeDMAT("../data/scene_data/glob_tree/rigs/skeleton_rig/surfaceW.DMAT", surface_W);
	Eigen::MatrixXd volume_W = surface_to_volume_weights(surface_W, bI, X, T);
	
	//igl::writeDMAT("../data/scene_data/glob_tree/rigs/skeleton_rig/surfaceW.DMAT", surface_W);
	//igl::writeDMAT("../data/scene_data/glob_tree/rigs/skeleton_rig/volumeW.DMAT", volume_W);
	Eigen::MatrixXd WT = volume_W.transpose();

	//finally
	Eigen::RowVectorXd colSums = volume_W.colwise().sum();
	Eigen::RowVectorXd rowSums = volume_W.rowwise().sum();
	this->W = volume_W;
	this->X = X;

	std::string volume_file_name = fs::path(surface_file_name).parent_path().string() + "/" + fs::path(surface_file_name).stem().string() + "_volume.json";
	printf("Conversion from surface skeleton rig to volume skeleton rig successful. Saving to %s ...\n", volume_file_name.c_str());
	write_rig_to_json(volume_file_name);

	init_rig_jacobian();
	init_null_space();
	S.resize(0, 3*X.rows());
	init_rig_selection_matrix(radius);
	//init_rig_selection_matrix(); // hold off on this for now
	//TODO: mesh bones into the mesh. hold off on this.


	//igl::writeDMAT(fs::path(surface_file_name).parent_path().string() + "/weights.dmat", W);
	//igl::writeDMAT(fs::path(surface_file_name).parent_path().string() + "/jacobian.dmat", J.toDense());
	//igl::writeDMAT(fs::path(surface_file_name).parent_path().string() + "/X.dmat", X);
	//igl::writeDMAT(fs::path(surface_file_name).parent_path().string() + "/T.dmat", T);
//	igl::writeMSH(fs::path(surface_file_name).parent_path().string() + "coarse_elephant.msh", X, T);
}

SkeletonRig::SkeletonRig(std::string rig_file_name, Eigen::MatrixXd& W, Eigen::MatrixXd& X, Eigen::MatrixXi& T, double radius)
{

	namespace fs = std::filesystem;

	json j;
	std::ifstream i(rig_file_name);
	i >> j;

	rig_pinning = j.value("rig_pinning", "ball");  //if htis is ball then we treat them like cylinders
	std::vector<std::vector<double>> X_list = j["vertices"];
	Eigen::MatrixXd surface_X;
	igl::list_to_matrix(X_list, surface_X); // hopefully this works
	Eigen::MatrixXd surface_W;
	read_bones_from_json(j, surface_X.rows(), surface_W, this->p0, this->pI, this->lengths);
	this->rig_type = "skeleton";
	this->W = W;
	this->X = X;
	init_rig_jacobian();
	init_rig_selection_matrix(radius);
}

SkeletonRig::SkeletonRig(std::string volume_file_name, double radius)
{
	read_rig_from_json(volume_file_name);

	init_rig_jacobian();
	init_rig_selection_matrix(radius);

}

bool SkeletonRig::read_rig_from_json(std::string filename)
{
namespace fs = std::filesystem;

	json j;

	std::ifstream i(filename);
	i >> j;

	std::string rig_pinning = j.value("rig_pinning", "ball");
	std::vector<std::vector<double>> X_list = j["X"];
	std::vector<std::vector<double>> W_list = j["W"];
	std::vector<double> p0_list = j["p0"];
	std::vector<double> length_list = j["lengths"];
	std::vector<int> pI_list = j["pI"];

	igl::list_to_matrix(X_list,X);
	igl::list_to_matrix(W_list,W);
	p0 = Eigen::Map<Eigen::VectorXd>(&p0_list[0], p0_list.size());
	pI = Eigen::Map<Eigen::VectorXi>(&pI_list[0], pI_list.size());
	lengths = Eigen::Map<Eigen::VectorXd>(&length_list[0], length_list.size());

	i.close();

	this->rig_type = "skeleton";
	return true;
}


bool SkeletonRig::write_rig_to_json(std::string filename)
{
    namespace fs = std::filesystem;

	filename = fs::path(filename).extension().string() == ".json" ? filename : filename + ".json";
	if (!fs::exists(fs::path(filename).parent_path()))
	{
		fs::create_directories(fs::path(filename).parent_path());
	}

	//write rig to json
	std::ofstream o(filename);
	json j;

	std::vector<std::vector<double>> X_list, W_list;
	std::vector<double> p0_list,  length_list;
	std::vector<int> pI_list;
        igl::matrix_to_list(X,X_list);
        igl::matrix_to_list(W,W_list);
	p0_list = std::vector<double>(p0.data(), p0.data() + p0.rows());
	pI_list = std::vector<int>(pI.data(), pI.data() + pI.rows());
	length_list = std::vector<double>(lengths.data(), lengths.data() + lengths.rows());

	j["rig_type"] = std::string("handle_rig");
	////save mesh positions
	j["X"] = X_list;
	j["W"] = W_list;
	j["p0"] = p0_list;
	j["pI"] = pI_list;
	j["lengths"] = length_list;
	o << j << std::endl;
	o.close();
	return true;
}
