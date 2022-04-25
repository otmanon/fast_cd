#include "rig.h"
#include "igl/cat.h"
#include "igl/repdiag.h"


AffineRig::AffineRig(Eigen::MatrixXd& x)
{
    X = x;
    dim = X.cols();
    n = X.rows();
    int vis_id = 1;
    gizmoPlugin = new igl::opengl::glfw::imgui::ImGuizmoWidget();  // bland gizmo pointer that will get destroyed as soon as we attach one
    init_parameters();
    init_jacobian();
    init_null_space();
    reset();
}

void AffineRig::init_parameters()
{
	Eigen::MatrixXd init;
	init.resize(3, 4);
	init.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
	init.col(3) = Eigen::Vector3d::Zero();
	//make sure this is row order flattened...nor sure if it is order flatten the first three rows
	p = Eigen::Map<Eigen::VectorXd>(init.data(), dim * (dim + 1), 1);
    p0 = p;
}

void AffineRig::init_jacobian()
{
	Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(X.rows(), 1);
	Eigen::MatrixXd X1 = igl::cat(2, X, ones);
	Eigen::MatrixXd J_dense;
	igl::repdiag(X1, dim, J_dense);
	J = (J_dense.sparseView());

}

void AffineRig::get_rig_motion(Eigen::VectorXd& p_curr)
{
	p_curr = p;
}

void AffineRig::init_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin)
{
    this->gizmoPlugin = gizmoPlugin;// remember pointer for later
    //fill out gizmo parameters. the gizmo itself is created and deleted outside of this object. It outlives this object
    //define rest transformation (centroid is actually at 0 displacement)
    rest_T = Eigen::MatrixXf::Identity(4, 4);
    //translate to center of mass
    rest_T.block(0, 3, 3, 1) = 0.5 * (X.colwise().maxCoeff() + X.colwise().minCoeff()).transpose().cast<float>();

    //set gizmos transformation to rest_T
    gizmoPlugin->T.block(0, 3, 3, 1) = rest_T.block(0, 3, 3, 1);
    gizmoPlugin->operation = ImGuizmo::TRANSLATE;
    gizmoPlugin->visible = true;
    //precalculatethe inverse
    rest_T_inv = rest_T.inverse();

    //unflattened rig parameters
    A = ((gizmoPlugin->T * rest_T_inv).block(0, 0, 3, 4).cast<double>()).transpose().eval();

    //starting rig parameters for an affine rig
    p = Eigen::Map<Eigen::VectorXd>(A.data(), 3 * 4, 1);    //flatten this
    gizmoPlugin->callback = [&](const Eigen::Matrix4f& T)
    {
        //get first 3 rows and 4 columns of gizmo transform
        A = ((T * rest_T_inv).block(0, 0, 3, 4).cast<double>()).transpose().eval();
        p = Eigen::Map<Eigen::VectorXd>(A.data(), 3 * 4, 1);    //flatten this tog et rig parameters
    };

}

void AffineRig::attach_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget* gizmoPlugin)
{
    //this->gizmoPlugin = *gizmoPlugin;
    //define rest transformation (centroid is actually at 0 displacement)
    rest_T = Eigen::MatrixXf::Identity(4, 4);
    //translate to center of mass
    rest_T.block(0, 3, 3, 1) = 0.5 * (X.colwise().maxCoeff() + X.colwise().minCoeff()).transpose().cast<float>();

    //set gizmos transformation to rest_T
    gizmoPlugin->T.block(0, 3, 3, 1) = rest_T.block(0, 3, 3, 1);
	gizmoPlugin->operation = ImGuizmo::TRANSLATE;
    //precalculatethe inverse
    rest_T_inv = rest_T.inverse();

    //unflattened rig parameters
    A = ((gizmoPlugin->T * rest_T_inv).block(0, 0, 3, 4).cast<double>()).transpose().eval();

    //starting rig parameters for an affine rig
    p = Eigen::Map<Eigen::VectorXd>(A.data(), 3 * 4, 1);    //flatten this
    gizmoPlugin->callback = [&](const Eigen::Matrix4f& T)
    {
        //get first 3 rows and 4 columns of gizmo transform
        A = ((T*rest_T_inv).block(0, 0, 3, 4).cast<double>()).transpose().eval();
        p = Eigen::Map<Eigen::VectorXd>(A.data(), 3 * 4, 1);    //flatten this tog et rig parameters
    };
}

void AffineRig::reset()
{
    //define rest transformation (centroid is actually at 0 displacement)
    Eigen::Matrix4f I = Eigen::MatrixXf::Identity(4, 4);
    rest_T = I;
    //translate to center of mass
    //replace to centroid of mesh, not C.O.M.
    rest_T.block(0, 3, 3, 1) = 0.5 * (X.colwise().maxCoeff() + X.colwise().minCoeff()).transpose().cast<float>();
    gizmoPlugin->T.block(0, 3, 3, 1) = rest_T.block(0, 3, 3, 1);

    Eigen::MatrixXd A = ((rest_T * rest_T.inverse()).block(0, 0, 3, 4).cast<double>()).transpose().eval();
    p = Eigen::Map<Eigen::VectorXd>(A.data(), 3 * 4, 1); // this one really hurts took way to long to find that this neede dto go here
    
}

bool AffineRig::key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier)
{
	switch (button)
	{
	 case 'r':
		 //change gizmo mode to rotations
		 gizmoPlugin->operation = ImGuizmo::ROTATE;
         return true;
	 case 'g': //g for grab like libigl. stay away from t as that'as used in default libigl viewer
	   //change gizmo mode to rotations
		 gizmoPlugin->operation = ImGuizmo::TRANSLATE;
		 return true;
     default:
         return false;
	}
}

void AffineRig::init_null_space()
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