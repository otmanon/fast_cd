#include "ConstraintControllerSquashStretch.h"
#include <igl/slice.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/centroid.h>
ConstraintControllerSquashStretch::ConstraintControllerSquashStretch(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo, float scale) :
	scale(scale), squash_factor(1)
{
    this->viewer = viewer;
    this->guizmo = guizmo;

    guizmo->visible = false;
    apply = true;
    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);

    //don't even care about surface I.
    //get left most indices
    double minX = X.col(0).minCoeff();
    double maxX = X.col(0).maxCoeff();
    double width = (maxX - minX);
    
    //get the 5% left most vertices
    double cutoff_min = minX + width * 0.1;
    double cutoff_max = maxX - width * 0.1;

    for (int i = 0; i < X.rows(); i++)
    {
        if (X(i, 0) < cutoff_min || X(i, 0) > cutoff_max)
        {
            bI.conservativeResize(bI.size() + 1);
            bI(bI.size() - 1) = i;
        }
    }
   // igl::unique(F, surfaceI);
    igl::slice(X, bI, 1, bc0);
    //it should only be the left most and right most ones
    bc = bc0;

    centroid = X.colwise().mean();
    //igl::centroid(X, T, centroid);
    
}

/*
Returns the desired boundary condition values of the vertices in question.
*/
void ConstraintControllerSquashStretch::get_interactive_motion(Eigen::MatrixXd& P)
{
    Eigen::MatrixXd rel_bc = bc0.rowwise() - centroid.transpose();

    Eigen::Matrix3d scale;
    scale << squash_factor, 0, 0,
        0, 1, 0,
        0, 0, 1;

    rel_bc = (scale * rel_bc.transpose()).transpose();
    bc = rel_bc.rowwise() + centroid.transpose();
    P = bc;
}

/*
Keycallback that controls the squash/stretch. Press A/a to squash and D/d to stretch
*/
bool ConstraintControllerSquashStretch::key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier)
{
    switch (button)
    {
    case 'A':
    case 'a':
        //why does this cause such a problem, whereas before it didnt?
       // data().set_face_based(!data().face_based);
        squash_factor += scale;
        return true;
    case 'D':
    case 'd':
        squash_factor -= scale;
        return true;
    }

}