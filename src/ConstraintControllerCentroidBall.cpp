#include "ConstraintControllerCentroidBall.h"
#include "igl/slice.h"
#include "igl/cat.h"
ConstraintControllerCentroidBall::ConstraintControllerCentroidBall(Eigen::MatrixXd& X, double r, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
	this->r = r;

	this->viewer = viewer;
	this->guizmo = guizmo;

	apply = true;
	guizmo->visible = true;
	
	Eigen::Vector3d center = X.colwise().mean();
	
	T = Eigen::Matrix4f::Identity();
	T.block(0, 3, 3, 1) = center.cast<float>();

	guizmo->T = T;
	guizmo->operation = ImGuizmo::TRANSLATE;
	T_rest_inv = T.inverse();

	//get vertices closest to central vertex
	for (int i = 0; i < X.rows(); i++)
	{
		const double d = (X.row(i).transpose() - center).squaredNorm();
		if (d < r * r)
		{
			bI.conservativeResize(bI.rows() + 1);
			bI(bI.rows() - 1) = i;
		}
	}
	igl::slice(X, bI, 1, bc0);

	this->guizmo->callback = [&](const Eigen::Matrix4f& A) {
		T = A * T_rest_inv;
	};
}

void ConstraintControllerCentroidBall::render(igl::opengl::glfw::Viewer& viewer)
{
	return;
}

void ConstraintControllerCentroidBall::get_interactive_motion(Eigen::MatrixXd& P)
{
	Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(bc0.rows(), 1);
	Eigen::MatrixXd V;
	igl::cat(2, bc0, ones, V);

	Eigen::MatrixXd P2;
	P2 = T.cast<double>() * V.transpose();
	P2.transposeInPlace();
	
	P = P2.block(0, 0, P2.rows(), 3);

}

