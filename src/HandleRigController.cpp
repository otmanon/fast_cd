#include "HandleRigController.h"
#include "igl/project.h"
#include "compute_handle_positions_from_rig_parameters.h"
#include "update_parameters_at_handle.h"
#include <igl/slice.h>
#include <igl/cat.h>

HandleRigMouseController::HandleRigMouseController(Eigen::VectorXd& p0,  igl::opengl::glfw::Viewer* viewer,igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
	// initialize guizmo to be at first rig parameters
	handleI = 0;
	this->p_rest = p0;
	this->p_rel = p0;  //these are relative parameters with respect to p_ref. For one rig, unflatten p_rest and p_rel to P_rest and P_rel. Then P_rel = P_rest.inverse() * P_deformed
	this->p_rest_inv = p0;

	vis_id = 1;
	Eigen::Matrix4f I, A_inv;
	I.setIdentity();
	for (int i = 0; i < p_rest.rows() / 12; i++)
	{ //update p_rel and p_rest in this loop.
		update_parameters_at_handle(this->p_rel, I, i);
		A_inv = matrix4f_from_parameters(p_rest, i);
		A_inv = A_inv.inverse().eval();
		update_parameters_at_handle(this->p_rest_inv, A_inv, i);
	}
	//get handle positiodsns
	compute_handle_positions_from_parameters(p_rest, V);

	this->guizmo = guizmo;
	init_guizmo_viewer(viewer, guizmo);
}

void HandleRigMouseController::draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{

	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		//	ImGui::Text("Rig Controller Menu");
			// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Rig Controller Menu", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		// Expose the same variable directly ...


		ImGui::End();
	};

}



void HandleRigMouseController::reset()
{
	Eigen::Matrix4f I, A_inv;
	I.setIdentity();
	for (int i = 0; i < p_rest.rows() / 12; i++)
	{ //update p_rel and p_rest in this loop.
		update_parameters_at_handle(this->p_rel, I, i);
		A_inv = matrix4f_from_parameters(p_rest, i);
		A_inv = A_inv.inverse().eval();
		update_parameters_at_handle(this->p_rest_inv, A_inv, i);
	}
	//get handle positiodsns
	compute_handle_positions_from_parameters(p_rest, V);

	this->guizmo->T = matrix4f_from_parameters(p_rest, handleI);

}
void  HandleRigMouseController::init_guizmo_viewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
//	this->guizmo = guizmo;
	vis_id = viewer->data_list.size();
	viewer->append_mesh();
	viewer->data_list[vis_id].clear();
	viewer->data_list[vis_id].set_points(V, Eigen::RowVector3d(0.8, 0.2, 0.2));
	viewer->data_list[vis_id].point_size = 50;


	this->guizmo->visible = true;
	this->guizmo->operation = ImGuizmo::TRANSLATE;
	this->guizmo->T = matrix4f_from_parameters(p_rest, handleI);


	this->guizmo->callback = [&](const Eigen::Matrix4f& A)
	{
		Eigen::Matrix4f A_inv = matrix4f_from_parameters(p_rest_inv, handleI);
		Eigen::Matrix4f T = A *A_inv;
		update_parameters_at_handle(p_rel, T, handleI);
	//	update_parameters_at_handle(p_rest, T, handleI);
		const int c = p_rest.rows() / 12;
		V.row(handleI) = A.block(0, 3, 3, 1).transpose().cast<double>();
		//compute_handle_positions_from_parameters(p, V);
		//V.row(handleI) = Eigen::RowVector3d(p(c*0 + handleI * 4), p(c*1 + handleI * 4), p(c * 2 + handleI * 4));
	};
}

void HandleRigMouseController::get_pinned_motion(Eigen::MatrixXd& X, std::vector<Eigen::VectorXi> bI, Eigen::VectorXd& bc)
{
	// used in satisfying S u = bc. bc must be ordered the same as the columns of S
	int ci = 0;
	Eigen::VectorXi I; 
	Eigen::MatrixXd BC_1_tmp, BC_tmp, BC;
	Eigen::MatrixXd X_b, X_b1;
	Eigen::Matrix4d T;

	Eigen::Matrix4f A_inv;
	Eigen::Matrix4f A;

	for (int b = 0; b < bI.size(); b++)
	{
		//get transofrmation matrix
		Eigen::Matrix4f A_inv = matrix4f_from_parameters(p_rest_inv, handleI);
		Eigen::Matrix4f A = matrix4f_from_parameters(p_rel, b);
		Eigen::Matrix4d T =( A ).cast<double>();

		I = bI[b];
		igl::slice(X, I, 1, X_b);

		Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(X_b.rows(), 1);
		igl::cat(2, X_b, ones, X_b1);

		BC_1_tmp = T * X_b1.transpose();
		BC_1_tmp.transposeInPlace();
		BC_tmp = BC_1_tmp.block(0, 0, BC_1_tmp.rows(), 3);

		BC.conservativeResize(BC.rows() + BC_tmp.rows(), 3);
		BC.bottomRows(BC_tmp.rows()) = BC_tmp - X_b;
	}

	bc = Eigen::Map<Eigen::VectorXd>(BC.data(), BC.rows()* 3);
}



void HandleRigMouseController::render(igl::opengl::glfw::Viewer& viewer)
{
	viewer.data_list[vis_id].set_points(V, Eigen::RowVector3d(0.8, 0.2, 0.2));
}



bool HandleRigMouseController::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	if (button == 2)
	{
		Eigen::Vector3d win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0.0);
		Eigen::Matrix4f view = viewer.core().view; Eigen::Matrix4f proj = viewer.core().proj; Eigen::Vector4f viewport = viewer.core().viewport;
		Eigen::MatrixXd p;
		igl::project(V, view, proj, viewport, p);
		Eigen::VectorXd D = (p.rowwise() - win.transpose()).rowwise().norm();
		int closest_tip;
		float min = D.minCoeff(&closest_tip);

		if (min < 100)
		{
			handleI = closest_tip;
			Eigen::Matrix4f T = Eigen::Matrix4f::Identity();
			
		//	for (int i = 0; i < T.rows()-1; i++)
			//	T.row(i) = p.middleRows(4 * handleI + 12*i, 4 * handleI + 4 + 12*handleI).cast<float>();
			Eigen::Matrix4f A = matrix4f_from_parameters(p_rel, handleI);
			Eigen::Matrix4f A_rest = matrix4f_from_parameters(p_rest, handleI);
			this->guizmo->T =  A_rest * A;
			return true;
		}
		else
		{
			return false;
		}
	}
	
	
	return false;
}

bool HandleRigMouseController::key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier)
{
	switch (button)
	{
	case 'g':
	case 'G':
		if (this->guizmo->operation == ImGuizmo::ROTATE)
			this->guizmo->operation = ImGuizmo::TRANSLATE;
		else
			this->guizmo->operation = ImGuizmo::ROTATE;
		return true;
		break;
	}
	return false;
}




Eigen::MatrixX4f HandleRigMouseController::matrix4f_from_parameters(Eigen::VectorXd& p, int ind)
{
	int c = 12;
	int num_p = p.rows() / c;

	Eigen::VectorXf a = p.middleRows(12 * ind, 12).cast<float>();
	Eigen::Matrix4f A = Eigen::Matrix4f::Identity();


	
	Eigen::MatrixXf T = Eigen::Map<Eigen::MatrixXf>(a.data(), 4, 3);
	A.topRows(3) = T.transpose();
	//for (int i = 0; i < 3; i++)
	//{
	//	Eigen::VectorXd row = p.middleRows(ind * 4 + i * 4 * (num_p),4);
	//	T.row(i) = row.transpose().cast<float>();
//
	//}

	return A;
}
