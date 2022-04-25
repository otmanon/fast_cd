#include "SkeletonRig.h"
#include "update_parameters_at_handle.h"
#include "surface_to_volume_weights.h"
#include "list_2D_to_matrix.h"
#include "matrix_to_2D_list.h"
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

#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <Eigen/Geometry>
std::string get_rig_file_format(std::string filename)
{
	//open surface_file as json
#ifdef WIN32
	namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif
	json j;
	std::ifstream i(filename);
	i >> j;

	if (j.value("format", "volume") == "volume")
	{
		return "volume";
	}
	else
	{
		return "surface";
	}

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
		bone_index = jbone["bone_idx"];  //maybe the bones are ordered in a different way in parent idx than what we loop throughh... so we need this
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
	

		if (wI.rows() == 0 || wV.rows() == 0)
		{
			std::cout << "Bone " + std::to_string(bone_index) + " has no effect on any vertices. This will result in singular rig Jacobian. Change your rig" << std::endl;
			exit(0);
		}
	}
}

SkeletonRig::SkeletonRig(std::string surface_file_name, Eigen::MatrixXd& X, Eigen::MatrixXi& T)
{

	printf("Converting a surface rig file format to a volume one. Will attempt to map from surface volume mesh X provided, to the coarse tet mesh we are trying to rig."
				"Then, we will diffuse surface weight values to volume weight values via a diffusion...\n");

		//open surface_file as json
	#ifdef WIN32
		namespace fs = std::filesystem;
	#else
		namespace fs = std::experimental::filesystem;
	#endif
	json j;
	std::ifstream i(surface_file_name);
	i >> j;


	rig_pinning = j.value("rig_pinning", "ball");  //if htis is ball then we treat them like cylinders
	std::vector<std::vector<double>> X_list = j["vertices"];
	Eigen::MatrixXd surface_X = list_2D_to_matrix(X_list); // hopefully this works
	Eigen::MatrixXd surface_W;
	read_bones_from_json(j, surface_X.rows(), surface_W, this->p0, this->pI, this->lengths);

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
		if (sums(i) > 1e-6)
		{
			sums(i) = 0;
			bI.conservativeResize(bI.rows() + 1);
			bI(bI.size() - 1) = i; //only get indices of non zero entries
		}
	}
	//Eigen::VectorXd W_tmp
	//igl::find(W_tmp, bI, bJ, bV);					//find all non-zero entries here, as they will now be treated as BOUNDARY conditions
	
	igl::slice(W_tmp, bI, 1, surface_W);
	Eigen::MatrixXd volume_W = surface_to_volume_weights(surface_W, bI, X, T);

	//finally
	this->W = volume_W;
	this->X = X;

	std::string volume_file_name = fs::path(surface_file_name).parent_path().string() + "/" + fs::path(surface_file_name).stem().string() + "_volume.json";
	printf("Conversion from surface skeleton rig to volume skeleton rig successful. Saving to %s ...\n", volume_file_name.c_str());
	write_rig_to_json(volume_file_name);

	init_rig_jacobian();
	init_null_space();
	S.resize(0, 3*X.rows());
	//init_rig_selection_matrix(); // hold off on this for now
	//TODO: mesh bones into the mesh. hold off on this.

}

SkeletonRig::SkeletonRig(std::string volume_file_name)
{
	read_rig_from_json(volume_file_name);

	init_rig_jacobian();
	init_rig_selection_matrix();

}

bool SkeletonRig::read_rig_from_json(std::string filename)
{
#ifdef WIN32
	namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif
	json j;

	std::ifstream i(filename);
	i >> j;

	rig_pinning = j.value("rig_pinning", "ball");
	std::vector<std::vector<double>> X_list = j["X"];
	std::vector<std::vector<double>> W_list = j["W"];
	std::vector<double> p0_list = j["p0"];
	std::vector<double> length_list = j["lengths"];
	std::vector<int> pI_list = j["pI"];

	X = list_2D_to_matrix(X_list);
	W = list_2D_to_matrix(W_list);
	p0 = Eigen::Map<Eigen::VectorXd>(&p0_list[0], p0_list.size());
	pI = Eigen::Map<Eigen::VectorXi>(&pI_list[0], pI_list.size());
	lengths = Eigen::Map<Eigen::VectorXd>(&length_list[0], length_list.size());
	return true;
}


bool SkeletonRig::write_rig_to_json(std::string filename)
{
#ifdef WIN32
	namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif
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
	X_list = matrix_to_2D_list(X);
	W_list = matrix_to_2D_list(W);
	p0_list = std::vector<double>(p0.data(), p0.data() + p0.rows());
	pI_list = std::vector<int>(pI.data(), pI.data() + pI.rows());
	length_list = std::vector<double>(lengths.data(), lengths.data() + lengths.rows());

	j["rig_type"] = "handle_rig";
	//save mesh positions
	j["X"] = X_list;
	j["W"] = W_list;
	j["p0"] = p0_list;

	o << std::setw(4) << j << std::endl;

	return true;
}
