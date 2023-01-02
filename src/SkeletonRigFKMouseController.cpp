#include "SkeletonRigFKMouseController.h"
#include "flattened_parameters_to_Affine3D_list.h"
#include "get_all_json_in_subdirs.h"
#include "load_animation.h"
#include <igl/forward_kinematics.h>

#include <igl/colon.h>
#include <igl/cat.h>
#include "update_parameters_at_handle.h"
#include "matrix4f_from_parameters.h"
#include <igl/project.h>

#include "get_tip_positions_from_parameters.h"
#include "get_joint_positions_from_parameters.h"
#include "get_relative_parameters.h"
#include "get_skeleton_mesh.h"


SkeletonRigFKMouseController::SkeletonRigFKMouseController(Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::VectorXd& pl, fast_cd_viewer* viewer, std::string animation_dir) : pl(pl), pI(pI)
{
		handleI = 0;
	thickness = 1e-2;
	this->p_rest = p0;
	int num_b = p_rest.rows() / 12;

	this->p_rel.resizeLike(p_rest);
	this->p_rest_inv.resizeLike(p_rest);
	//convert flattened matrix to affine matrix
	//flattened_parameters_to_Affine3d_list(this->p0, A0);

	//Get edge connectivity  YIKES can't believe I wrote this...How did this work???
	Eigen::VectorXi I1, I2;
	igl::colon(0, num_b - 1, I1);
	I2 = I1.array() + num_b;
	BE.resize(num_b, 2);
	BE.col(0) = I1;
	BE.col(1) = I2;

	reset();

	init_guizmo_viewer(viewer->igl_v, viewer->guizmo);


	std::string custom_anim_name = "custom_anim";
	recording= false;
	current_animation_id = -1;
	loaded_anim = false;
	pause = true;

	//this->animation_dir = animation_dir;
	//get_all_json_in_dir(animation_dir, animation_filepaths, animation_filenames);
	//anim_step = 0;
}

SkeletonRigFKMouseController::SkeletonRigFKMouseController(Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::VectorXd& pl, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo, std::string animation_dir)
	: pl(pl), pI(pI)
{
	handleI = 0;
	thickness = 1e-2;
	this->p_rest = p0;
	int num_b = p_rest.rows() / 12;

	this->p_rel.resizeLike(p_rest);
	this->p_rest_inv.resizeLike(p_rest);
	//convert flattened matrix to affine matrix
	//flattened_parameters_to_Affine3d_list(this->p0, A0);

	//Get edge connectivity  YIKES can't believe I wrote this...How did this work???
	Eigen::VectorXi I1, I2;
	igl::colon(0, num_b - 1, I1);
	I2 = I1.array() + num_b;
	BE.resize(num_b, 2);
	BE.col(0) = I1;
	BE.col(1) = I2;

	reset();

	init_guizmo_viewer(viewer, guizmo);


	std::string custom_anim_name = "custom_anim";
	recording= false;
	current_animation_id = -1;
	loaded_anim = false;
	pause = true;
	this->animation_dir = animation_dir;
	get_all_json_in_dir(animation_dir, animation_filepaths, animation_filenames);
	anim_step = 0;

	
}


//show the current skeleton visually
void SkeletonRigFKMouseController::init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
	this->guizmo = guizmo;
	vis_id = 2; // all skeleton rigs by default go to 2 here
	while (viewer->data_list.size() <= 2)
		viewer->append_mesh();
	viewer->data_list[vis_id].clear();
	
	//bone radius


	viewer->data_list[vis_id].point_size = 10;
	viewer->data_list[vis_id].line_width = 10;

	this->guizmo->visible = true;
	this->guizmo->operation = ImGuizmo::ROTATE;

	//when we pick the handle, make sure it is aligned to that bones rotation matrix, instead of just the identity
	this->guizmo->T = matrix4f_from_parameters(p_rest, handleI);


	this->guizmo->callback = [&](const Eigen::Matrix4f& A)
	{
		int num_b = p_rest.rows() / 12;
		Eigen::Matrix4f A_inv = matrix4f_from_parameters(p_rest_inv, handleI);
		Eigen::Matrix4f T = (A * A_inv);  // P_rel... but should be P_rel relative to parent

		Eigen::Matrix3d R = T.block(0, 0, 3, 3).cast<double>();
		R_rel_parent[handleI] = Eigen::Quaterniond(R);
		Eigen::MatrixXd A_list;

		//need to compute crucial R_rel_parent at current bone only.


		igl::forward_kinematics(C0, BE, pI, this->R_rel_parent, A_list); // A_list contains global transformation matrices for each bone
		
	    //From A_list, need to flatten it to sol_p.
		Eigen::VectorXd p;
		p.resize(12 * num_b);
		Eigen::Matrix4f AG = Eigen::Matrix4f::Identity();
		for (int i = 0; i < num_b; i++)
		{
			AG.block(0, 0, 3, 4) = A_list.block(i * 4, 0, 4, 3).cast<float>().transpose();
			update_parameters_at_handle(p, AG, i);
		}

		p_rel = p;
		Eigen::VectorXd p_glob;
		//from sol_p, get relative rig parameters... this is what CD ultimately uses!
	//	get_relative_parameters(p_rest, sol_p, p_rel);
		get_absolute_parameters(p_rest, p, p_glob);
		//get new C... for visualization purposes. 
		Eigen::MatrixXd joints, tips;
		get_tip_positions_from_parameters(p_glob, this->pl, tips);
		get_joint_positions_from_parameters(p_glob, joints);
		C.resize(joints.rows() + tips.rows(), 3);
		C.topRows(joints.rows()) = joints;
		C.bottomRows(tips.rows()) = tips;

		get_skeleton_mesh(thickness, p_glob, this->pl, renderV, renderF, renderC);
		//Just nee
	//	get_tip_positions_from_parameters(p_rel)
		//V.row(handleI) = A.block(0, 3, 3, 1).transpose().cast<double>();
		//
		//compute_handle_positions_from_parameters(sol_p, V);
		//V.row(handleI) = Eigen::RowVector3d(sol_p(c*0 + handleI * 4), sol_p(c*1 + handleI * 4), sol_p(c * 2 + handleI * 4));
	};
}




void SkeletonRigFKMouseController::draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{
	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
	//	ImGui::Text("Rig Controller Menu");
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Rig Controller Menu", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		ImGui::SliderFloat("Bone Thickness", &thickness, 1e-3, 1);


		if (ImGui::CollapsingHeader("Animation"))
		{
			//List all the animations in our directory
			if (ImGui::BeginListBox("Animations"))
			{
				for (int n = 0; n < animation_filenames.size(); n++)
				{
					const bool is_selected = (current_animation_id == n);
					if (ImGui::Selectable(animation_filenames[n].c_str(), is_selected))
					{
						//start animation!!!
						current_animation_id = n;
						anim_step = 0;
						load_animation_and_fit_to_rest_skeleton_json(animation_filepaths[n], p_rest, this->pI,  anim_P, is_global_anim);
						loaded_anim = true;
						guizmo->visible = false;
					
					}
						//new_animation = n;

					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndListBox();
			}

			if (ImGui::Button("Pause/Play Anim"))
			{
				pause = !pause;
			}

			if (ImGui::Button("Show/Hide Guizmo"))
			{
				guizmo->visible = !guizmo->visible;
			}
		}

		if (ImGui::CollapsingHeader("Record Rig"))
		{
			if (ImGui::Button("Pause/Play Anim"))
			{
				pause = !pause;
			}



		}
		ImGui::End();
	};

}


//show your mesh girl!
void SkeletonRigFKMouseController::render(igl::opengl::glfw::Viewer& viewer)
{
	viewer.data_list[vis_id].set_mesh(renderV, renderF);	
	viewer.data_list[vis_id].set_colors(renderC);
	viewer.data_list[vis_id].set_edges(C, BE, Eigen::RowVector3d(0, 0, 0));
	viewer.data_list[vis_id].set_points(C, Eigen::RowVector3d(0.8, 0.2, 0.2));
	viewer.data_list[vis_id].point_size = 5;
	viewer.data_list[vis_id].line_width = 5;



}




void SkeletonRigFKMouseController::reset()
{
	int num_b = p_rest.rows() / 12;
	//get skeleton geometry (concatenate joints with tips)
	Eigen::MatrixXd joints, tips;
	get_tip_positions_from_parameters(this->p_rest, this->pl, tips);
	get_joint_positions_from_parameters(this->p_rest, joints);
	igl::cat(1, joints, tips, C);
	C0 = C;
	Eigen::Matrix4f I, A_inv, A_rest;
	I.setIdentity();
	for (int i = 0; i < num_b; i++)
	{ //update p_rel and p_rest in this loop.
		update_parameters_at_handle(this->p_rel, I, i);
		A_rest = matrix4f_from_parameters(p_rest, i);
		A_inv = A_rest.inverse();
		update_parameters_at_handle(this->p_rest_inv, A_inv, i);
	}
	Eigen::VectorXd p_glob;
	get_absolute_parameters(p_rest, p_rel, p_glob);
	get_skeleton_mesh(thickness, p_glob, this->pl, renderV, renderF, renderC);

	//get our local rotations for each bone...
	Eigen::Quaterniond q = Eigen::Quaterniond::Identity();
	R_rel_parent.resize(num_b, q);

}



void SkeletonRigFKMouseController::set_scripted_motion(int step) 
{
	if (loaded_anim && !pause)
	{
		
		Eigen::VectorXd p_tmp = anim_P.col(anim_step);
		get_relative_parameters(p_rest, p_tmp, p_rel);
		get_skeleton_mesh(thickness, p_tmp, this->pl, renderV, renderF, renderC);

		Eigen::MatrixXd joints, tips;
		get_tip_positions_from_parameters(p_tmp, this->pl, tips);
		get_joint_positions_from_parameters(p_tmp, joints);
		C.resize(joints.rows() + tips.rows(), 3);
		C.topRows(joints.rows()) = joints;
		C.bottomRows(tips.rows()) = tips;
		anim_step += 1;
		anim_step = anim_step % anim_P.cols();
	}


};

bool SkeletonRigFKMouseController::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	if (button == 2)
	{
		Eigen::Vector3d win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0.0);
		Eigen::Matrix4f view = viewer.core().view; Eigen::Matrix4f proj = viewer.core().proj; Eigen::Vector4f viewport = viewer.core().viewport;
		Eigen::MatrixXd p;

		Eigen::MatrixXd joints; Eigen::VectorXd p_glob;
		get_absolute_parameters(p_rest, p_rel, p_glob);
		get_joint_positions_from_parameters(p_glob, joints);
		igl::project(joints, view, proj, viewport, p);
		Eigen::VectorXd D = (p.rowwise() - win.transpose()).rowwise().norm();
		int closest_tip;
		float min = D.minCoeff(&closest_tip);

		if (min < 100)
		{
			handleI = closest_tip;
			Eigen::Matrix4f T = Eigen::Matrix4f::Identity();

			//	for (int i = 0; i < T.rows()-1; i++)
				//	T.row(i) = sol_p.middleRows(4 * handleI + 12*i, 4 * handleI + 4 + 12*handleI).cast<float>();
			Eigen::Matrix4f A = matrix4f_from_parameters(p_glob, handleI);
			this->guizmo->T = A;
			return true;
		}
		else
		{
			return false;
		}
	}


	return false;
}