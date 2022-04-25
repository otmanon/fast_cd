#include "lbs_rig.h"
#include <igl/colon.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include "projected_gradient_descent.h"
#include "levenberg_marquardt.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include "interweaving_matrix.h"
#include <igl/centroid.h>

#include "Bone.h"
#include "euler_angles_to_transform.h"
#include "forward_kinematics.h"
#include "copy_skeleton_at.h"
#include "transformed_tips.h"
#include "kinematics_jacobian.h"
#include "end_effectors_objective_and_gradient.h"
#include "read_anim_file.h"

#include "skeleton_visualization_mesh.h"
#include "get_all_json_in_subdirs.h"

#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
LBSRig::LBSRig(std::string filename, Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::VectorXi& bI , std::string anim_file_dir)
{
#ifdef WIN32
	namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif
	//classic
	this->X = X;
	this->bI = bI;

	//rig is always set to 1. simulation hook is set to 0
	vis_id = 1;
	//open our json file!
	std::ifstream infile(filename);
	if (!infile)
	{
		std::cout << "No rig json file found at " + filename << std::endl;
		return;
	}
	json j;
	infile >> j;

	surfaceX = read_vertices_from_json(j);
	
	surfaceF = read_faces_from_json(j);

	skeleton = read_bones_from_json(j);


	fit_bone_transforms_to_input_X();

	
	thickness = 0.01/best_fit_scale.block(0, 0, 3, 3).diagonal().mean();
	double maxV = X.col(1).maxCoeff();
	double minV = X.col(1).minCoeff();
	//should assert here that surfaceX is indeed X[bI] == surfaceX
	Eigen::MatrixXd surfaceXtest;
	igl::slice(X, bI, 1, surfaceXtest);
	double length_scale = (surfaceX.rowwise().maxCoeff() - surfaceX.rowwise().minCoeff()).mean();
	Eigen::MatrixXd error = (surfaceX - surfaceXtest).cwiseAbs()/length_scale; //this should be relative to the size of the mesh
	assert(error.maxCoeff() < 1e-5 && "surface vertices in json file do \
		not match with X[bI]. If using blender, ensure when importing .obj \
		file that you retain the vertex order, and that the mesh you give it\
		 satisfies meshX = volumeX[bI]. \n bI here is obtained through \
		igl::unique(igl::boundary_facets(T)). It  seems the surface mesh\
		 provided by tetWild automatically satisfies this.");


	init_parameters();


	surfaceW = weight_matrix_from_bones(skeleton);

	Eigen::VectorXd row_sum = surfaceW.rowwise().sum();
	Eigen::VectorXd col_sum = surfaceW.colwise().sum();
	W = surface_to_volume_weights(surfaceW, bI, X, T);

	init_jacobian();

	P = get_all_bone_heads(skeleton);

	mouse_dragging = false;
	added_handle = false;
	moved_handle = false;


	interaction_mode = INTERACTION_MODE::INVERSE_KINEMATICS;
	ik_solver = ik::IK_SOLVER::LM;
	max_iter = 10;
	mu = 1e-5;
	tol_exp = 10;
	tol = pow(10, -tol_exp);

	INTERACTION_MODE new_interaction_mode;
	int frame_i = 0;
	int total_frames = anim_T.size();

	skeleton_visualization_mesh(skeleton, thickness, renderV, renderF, renderC);

	//
	anim_paths.clear(); anim_names.clear();
	get_all_json_in_dir(fs::path(filename).parent_path().string() + "/anim/", anim_paths, anim_names);
	current_anim_id = -1;
	new_anim_id = -1;
}

void LBSRig::reset()
{
	handleI.resize(0);
	handleV.resize(0, 3);
	b.resize(0);
	xb0.resize(0);


	for (int i = 0; i < skeleton.size(); i++)
	{
		skeleton[i].xzx = Eigen::Vector3d(0, 0, 0);
	}
	frame_i = 0;

	p = p0;
} 

void LBSRig::init_jacobian()
{
	double w_ij;
	Eigen::Vector3d vi, wvi;
	int b;
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(12 * X.rows() * num_b); // 12 entries for each bone-vertex pair
	int v_dim = 3;
	int a_dim = 4;
	int row, col;
	
	for (int i = 0; i < X.rows(); i++)
	{
		vi = X.row(i);
		for (int j = 0; j < num_b; j++)
		{
			b = skeleton[j].weight_index;
			w_ij = W(i, b);
			wvi = w_ij*vi;
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
	J.resize(v_dim * X.rows(), num_p);
	J.setFromTriplets(tripletList.begin(), tripletList.end());

	Eigen::SparseMatrix<double> S;
	interweaving_matrix(X.rows(), X.cols(), S);
	J = S * J;		//output is row order flattened, we want column order flattened
}

void LBSRig::init_parameters() 
{
	num_b = skeleton.size();
	num_p = num_b * 12;
	std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > T;
	forward_kinematics(skeleton, T); //may need to multiply each T by the rest_bones.
	
	//ALSO make sure each bone.weight_index is equivalent to its place in the bone_list
	update_rig_parameters(T);
	p0 = p;
};

//TODO think a bit harder about the arguments being passed here... should we be storing this T as a global list, Everyone uses it!
void LBSRig::update_rig_parameters(std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>& T)
{
	p.conservativeResize(num_p);

	int rows = 3;
	int cols = 4;
	Eigen::Affine3d T_total;
	//do row wise flattening of rig parameters
	for (int i = 0; i < T.size(); i++)
	{
		T_total = T[i]; // *skeleton[i].rest_T;
		for (int row = 0; row < T[i].rows(); row++)
		{
			p.block(i * cols * rows + row * cols, 0, cols, 1) = T_total.matrix().row(row).transpose();
		}
	}
}

Eigen::MatrixXd LBSRig::surface_to_volume_weights(Eigen::MatrixXd& surfaceW, Eigen::VectorXi& bI, Eigen::MatrixXd& X, Eigen::MatrixXi& T)
{
	Eigen::MatrixXd W;
	Eigen::SparseMatrix<double> C, Aeq;
	igl::cotmatrix(X, T, C);
	C *= -1;
	Eigen::MatrixXd Beq;

	Eigen::MatrixXi BI = bI.replicate(1, surfaceW.rows());
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(C.rows(), surfaceW.cols());
	igl::min_quad_with_fixed_data<double> data;
	igl::min_quad_with_fixed_precompute(C, bI, Aeq, true, data);
	igl::min_quad_with_fixed_solve(data, B, surfaceW, Beq, W);
	//normalize weight matrix
	Eigen::VectorXd sums = W.rowwise().sum();
	W.array().colwise() /= W.rowwise().sum().array();
	return W;
}

Eigen::MatrixXd LBSRig::read_vertices_from_json(json& j)
{
	//get vertices. Ensure that the vertices are the first #vertices of the simulation V.
	std::vector<std::vector<double>> vertex_list = j["vertices"];
	Eigen::MatrixXd V(vertex_list.size(), 3);
	for (int i = 0; i < vertex_list.size(); i++)
		V.row(i) = Eigen::Map<Eigen::RowVector3d>(vertex_list[i].data(), 3);
	return V;
}

Eigen::MatrixXi LBSRig::read_faces_from_json(json& j)
{
	//get faces. 
	std::vector<std::vector<int>> face_list = j["faces"];
	Eigen::MatrixXi F(face_list.size(), 3);
	for (int i = 0; i < face_list.size(); i++)
		F.row(i) = Eigen::Map<Eigen::RowVector3i>(face_list[i].data(), 3);

	return F;
}

Skeleton LBSRig::read_bones_from_json(json& j)
{
	const int num_bones = j["bones"].size();
	Skeleton bone_list;
	bone_list.reserve(num_bones);
	//use this to convert from flattened matrix to 4x4 affine matrix
	auto parse_affine = [](const json& j) -> Eigen::Affine3d
	{
		return Eigen::Affine3d((Eigen::Matrix4d(4, 4) <<
			j[0][0], j[0][1], j[0][2], j[0][3],
			j[1][0], j[1][1], j[1][2], j[1][3],
			j[2][0], j[2][1], j[2][2], j[2][3],
			0, 0, 0, 1).finished());
	};

	// parse bones to initialize skeleton
	int parent_index, bone_index;
	double length;
	Eigen::VectorXi wI;
	Eigen::VectorXd wV;
	for (const json& jbone : j["bones"])
	{
		parent_index = jbone["parent_idx"];
		bone_index = jbone["bone_idx"];
		Eigen::Affine3d A = parse_affine(jbone["bone_transform"]);
		double length = jbone["bone_length"];

		std::vector<int> w_i = jbone["vert_indeces_of_bone"];
		wI = Eigen::Map<Eigen::VectorXi>(w_i.data(), w_i.size());
		std::vector<double> w_v = jbone["weights_of_bone"];
		wV = Eigen::Map<Eigen::VectorXd>(w_v.data(), w_v.size());

		bone_list.emplace_back(
			parent_index,
			bone_index,
			A,
			length, wI, wV);

		if (wI.rows() == 0 || wV.rows() == 0)
		{
			std::cout << "Bone " + std::to_string(bone_index) + " has no effect on any vertices. This will result in singular rig Jacobian. Change your rig" << std::endl;
			exit(0);
		}
	}

	return bone_list;
}

Eigen::MatrixXd LBSRig::weight_matrix_from_bones(Skeleton& bone_list)
{
	// initialize weight matrix
	Eigen::MatrixXd W = Eigen::MatrixXd::Zero(surfaceX.rows(), bone_list.size());

	int i = 0;
	Eigen::VectorXi wI;
	Eigen::VectorXd wV;
	for (const Bone& bone : bone_list)
	{
		wI = bone.wI;
		wV = bone.wV;
		for (int index = 0; index < wI.rows(); index++)
		{
			W(wI[index], bone.weight_index) = wV(index);
		}
	}

	//normalize weight matrix
	Eigen::VectorXd sums = W.rowwise().sum();
	W.array().colwise() /= W.rowwise().sum().array();

	Eigen::VectorXd bone_cont = W.colwise().sum();
	return W;
}

Eigen::MatrixXd LBSRig::get_all_bone_heads(Skeleton& bone_list)
{
	Eigen::VectorXi bI;
	igl::colon(0, num_b-1, bI);		//hate that the stop value is inclusive
	Eigen::VectorXd P_vec = transformed_tips(bone_list, bI);
	Eigen::MatrixXd P = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(P_vec.data(), P_vec.rows() / 3, 3);
	return P;
}

//called in simulation thread... don't like that we are doing inverse kinematics and everything here. TODO: clean this method up, get_rig_motion should only return the rig parameters. 
//Other fancy lbs stuff can be done somewhere else. 
void LBSRig::get_rig_motion(Eigen::VectorXd& p_curr) 
{
	poll_rig_changes();

	//these are our rig parameters... we recalculate them each timestep
	std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > T;
	Eigen::VectorXd x;
	x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
	switch (interaction_mode)
	{
	case INTERACTION_MODE::INVERSE_KINEMATICS:
		if (b.rows() > 0)
		{
			inverse_kinematics(skeleton, b, xb0, max_iter, tol, ik_solver, mu);
		}

		forward_kinematics(skeleton, T);

		update_rig_parameters(T);

		p_curr = this->p;
		//keep the head tips updates for the viewer'as sake
		P = get_all_bone_heads(skeleton);
		skeleton_visualization_mesh(skeleton, thickness, renderV, renderF, renderC);

		break;
	case INTERACTION_MODE::ANIMATION:

		if (current_anim_id > -1)
		{
		

			int i = frame_i / 2;
			std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> frame_T = anim_T[i];

			update_rig_parameters(frame_T);
			p_curr = p;

			frame_i += 1;
			frame_i = frame_i % (2 * anim_T.size());

			skeleton_visualization_mesh(skeleton, frame_T, thickness, renderV, renderF, renderC);
		}
		else {
			p_curr.resize(J.rows());
			p_curr.setZero();
			std::cout << "Didn't select a rig animation! Click on a rig animation in the animation list!" << std::endl;
		}

		break;
	}
	
	
};


void LBSRig::fit_bone_transforms_to_input_X()
{


	Eigen::RowVector3d s = (X.colwise().maxCoeff() - X.colwise().minCoeff()).array() / (surfaceX.colwise().maxCoeff() - surfaceX.colwise().minCoeff()).array();
	Eigen::Matrix4d S;
	S << s(0), 0, 0, 0,
		0, s(1), 0, 0,
		0, 0, s(2), 0,
		0, 0, 0, 1;

	best_fit_scale = S;
	//change rest T of each bone
	Eigen::Affine3d T;
	for (int i = 0; i < skeleton.size(); i++)
	{
		//we extract the transform
		T = skeleton[i].rest_T;
		//we fix it with ourscaling matrix
		T.matrix() = S * T.matrix();
		//we place the transform back into its bone
		skeleton[i].rest_T = T; //like surgery :-) 
	}
	//update surfaceX so that it'as corrext scale
	surfaceX = (S.block(0, 0, 3, 3) * surfaceX.transpose().eval()).transpose().eval();
	
	//now that we have scale let'as fit translation!
	Eigen::MatrixXd surfaceXtest;
	igl::slice(X, bI, 1, surfaceXtest);

	Eigen::RowVector3d surfaceG = surfaceX.colwise().mean();
	Eigen::RowVector3d realG = surfaceXtest.colwise().mean();
	Eigen::RowVector3d diffG = realG - surfaceG;

	best_fit_translation = diffG;
	for (int i = 0; i < skeleton.size(); i++)
	{
		skeleton[i].rest_T.matrix().block(0, 3, 3, 1) += diffG.transpose();
	}

	surfaceX.rowwise() += diffG;


}

void LBSRig::fit_anim_transforms_to_input_X()
{
	for (int fi = 0; fi < anim_T.size(); fi++)
	{
		for (int bi = 0; bi < anim_T[fi].size(); bi++)
		{
			//we extract the transform
			
			//we fix it with ourscaling matrix
			anim_T[fi][bi].matrix() = best_fit_scale * anim_T[fi][bi].matrix();

			anim_T[fi][bi].matrix().block(0, 3, 3, 1) += best_fit_translation;
		}
	}
}


void LBSRig::init_viewer(igl::opengl::glfw::Viewer& viewer)
{
	this->viewer = &viewer;
	viewer.append_mesh();
	viewer.data_list[vis_id].clear();
	viewer.data_list[vis_id].set_mesh(renderV, renderF);
	viewer.data_list[vis_id].set_face_based(true);
	viewer.data_list[vis_id].set_colors(renderC);

	viewer.data_list[vis_id].point_size = 50;
	//viewer.data_list[vis_id].add_points(P, blue);
}

void LBSRig::render(igl::opengl::glfw::Viewer& viewer)
{
	//reload visualization render data
	
	viewer.data_list[vis_id].set_mesh(renderV, renderF);
	viewer.data_list[vis_id].set_face_based(true);
	viewer.data_list[vis_id].set_colors(renderC);

	if (handleI.rows() > 0)
	{
		viewer.data_list[vis_id].set_points(handleV, handleC);
	}
	else
		viewer.data_list[vis_id].clear_points();

}

void LBSRig::poll_rig_changes()
{
	if (added_handle)
	{
		//for some reason map doesn't like to be given a transpose()
		Eigen::MatrixXd handleV2 = handleV.transpose();
		xb0 = Eigen::Map<Eigen::VectorXd>(handleV2.data(), handleV.rows() * 3);
		b = handleI;

		added_handle = false;
	}
	if (moved_handle)
	{
		//for some reason map doesn't like to be given a transpose()
		Eigen::MatrixXd handleV2 = handleV.transpose();
		xb0 = Eigen::Map<Eigen::VectorXd>(handleV2.data(), handleV.rows() * 3);
		moved_handle = false;
	}
	if (current_anim_id != new_anim_id)
	{
		std::cout << "reading rig animation " << anim_names[new_anim_id] << "..." << std::endl;
		std::string anim_file_dir = anim_paths[new_anim_id];
		anim_T = read_anim_file(anim_file_dir, skeleton, best_fit_scale, best_fit_translation);
		fit_anim_transforms_to_input_X();
		update_rig_parameters(anim_T[0]);
		current_anim_id = new_anim_id;
		frame_i = 0;

	}
}

void LBSRig::draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{
	if (ImGui::CollapsingHeader("Rig Menu"))
	{
		ImGui::Text("Rig Mode");
		if (ImGui::RadioButton("Animation", (int*)&interaction_mode, (int)INTERACTION_MODE::ANIMATION))
		{
			reset();
		}
		if (ImGui::RadioButton("IK", (int*)&interaction_mode, (int)INTERACTION_MODE::INVERSE_KINEMATICS))
		{
			reset();
		}
		if (interaction_mode == INTERACTION_MODE::INVERSE_KINEMATICS)
		{
			if (ImGui::CollapsingHeader("Inverse Kinematics"))
			{
				ImGui::Text(" IK Solver");
				ImGui::RadioButton("Gradient Descent", (int*)&ik_solver, (int)ik::IK_SOLVER::GD);
				ImGui::RadioButton("Levenberg Marquardt", (int*)&ik_solver, (int)ik::IK_SOLVER::LM);
				if (ik_solver == ik::IK_SOLVER::LM)
				{
					ImGui::SliderFloat("mu", (float *) &mu, 1e-12, 1e2, "%.3g", ImGuiSliderFlags_Logarithmic);
				}
				ImGui::InputInt("max iter", &max_iter, 10, 50);
				ImGui::Text("tol: 1e-"); ImGui::SameLine();
				ImGui::SliderInt("", &tol_exp, 1, 16); ImGui::SameLine();
				tol = pow(10, -tol_exp);
			}

		}
		
		if (ImGui::CollapsingHeader("Anims"))
		{
			if (ImGui::BeginListBox("Anims"))
			{
				for (int n = 0; n < anim_names.size(); n++)
				{
					const bool is_selected = (current_anim_id == n);
					if (ImGui::Selectable(anim_names[n].c_str(), is_selected))
						new_anim_id = n;

					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndListBox();
			}
		}
		
	}
}

bool LBSRig::mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	if (mouse_dragging)
	{
		handleC.row(currentHandle) = blue;
		currentHandle = -1;
		viewer.data_list[vis_id].set_points(handleV, handleC);
	}
	mouse_dragging = false;


	return false;
}

bool LBSRig::mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	if (mouse_dragging)
	{
		mouse_drag_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, mouse_win(2));
		mouse_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, mouse_win(2));

		igl::unproject(
			mouse_drag_win,
			viewer.core().view,
			viewer.core().proj,
			viewer.core().viewport,
			mouse_drag_world);


		igl::unproject(
			mouse_win,
			viewer.core().view,
			viewer.core().proj,
			viewer.core().viewport,
			mouse_world);


		handleV.row(currentHandle) = mouse_world;
		handleC.row(currentHandle) = red;


		

		viewer.data_list[vis_id].set_points(handleV, handleC);
		moved_handle = true;
		return true;
	}
	else
	{
		return false;
	}
}

bool LBSRig::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	if (button == 2) //right click means we are making a new vertex
	{
		Eigen::Vector2f mouse_win_2d = Eigen::Vector2f(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y);
		mouse_win = Eigen::Vector3d(mouse_win_2d.x(), mouse_win_2d.y(), 0);
		
		Eigen::MatrixXd CP; //closest point
		igl::project(
			P,
			viewer.core().view,
			viewer.core().proj, viewer.core().viewport, CP);
		Eigen::VectorXd D = (CP.rowwise() - mouse_win.transpose()).rowwise().norm();
		int closest_tip;
		float min = D.minCoeff(&closest_tip);
		if (closest_tip != -1)
		{
			for (int i = 0; i < handleI.rows(); i++)
			{
				if (handleI(i) == closest_tip) return false;
			}
			Eigen::RowVector3d newV = P.row(closest_tip);
		//	viewer.data_list[vis_id].add_points(newV, red);
			handleV.conservativeResize(handleV.rows() + 1, 3);
			handleV.row(handleV.rows() - 1) = newV;
			handleC.conservativeResize(handleC.rows() + 1, 3);
			handleC.row(handleC.rows() - 1) = blue;
			handleI.conservativeResize(handleI.rows() + 1);
			handleI(handleI.rows() - 1) = closest_tip;

			viewer.data_list[vis_id].set_points(handleV, handleC);
			added_handle = true;
			return true;
		}
	}
	if (button == 0) //left click means we are moving a vertex
	{   //copying geometry processing assignment, 
		if (handleI.rows() > 0)
		{
			mouse_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
			Eigen::MatrixXd CP; //closest point
			igl::project(
				handleV,
				viewer.core().view,
				viewer.core().proj, viewer.core().viewport, CP);
			Eigen::VectorXd D = (CP.rowwise() - mouse_win.transpose()).rowwise().norm();
			currentHandle = (D.minCoeff(&currentHandle) < 50) ? currentHandle : -1;
			if (currentHandle != -1)
			{
				mouse_dragging = true;
				mouse_win(2) = CP(currentHandle, 2);
				return true;
		
			}
		}
	}
	return false;
}


