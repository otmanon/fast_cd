#pragma once
#include "fast_cd_viewer.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
struct fast_cd_viewer_vertex_selector : public fast_cd_viewer {

	int currI;
	MatrixXd C;
	VectorXi CI;

    int vid; // which of the data layers should the viewer be configured to
    RowVector3f last_mouse;

    //these aren't needed for the visuals to work, but are useful for people to detect what happend
    //this timestep

    bool vertex_added;
	fast_cd_viewer_vertex_selector(int vid = 0) : fast_cd_viewer()
	{
        this->vid = vid;
        vertex_added = false;
        ////////////// MOUSE CLICK INTERACTION ////////////////
        currI = -1;  
        last_mouse = RowVector3f::Zero();
        igl_v->callback_mouse_down =
            [&](igl::opengl::glfw::Viewer&, int button, int action)->bool
        { 
            last_mouse = Eigen::RowVector3f(
                igl_v->current_mouse_x, igl_v->core().viewport(3) - igl_v->current_mouse_y, 0);
            if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Left)
                printf("Left Click\n");
            
            else if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Right)
                printf("Right Click \n");//
            
            printf("button %d \n", button);
            printf("action %d \n", action);
            if (igl::opengl::glfw::Viewer::MouseButton(button) == igl::opengl::glfw::Viewer::MouseButton::Right)  //if right click, we are placing handles
            {
                // Find closest point on mesh to mouse position
                int fid;
                Eigen::Vector3f bary;
                MatrixXd U = igl_v->data_list[this->vid].V;
                MatrixXi F = igl_v->data_list[this->vid].F;
                if (igl::unproject_onto_mesh(
                    last_mouse.head(2),
                    igl_v->core().view,
                    igl_v->core().proj,
                    igl_v->core().viewport,
                    U, F,
                    fid, bary))
                {
                    long c;
                    bary.maxCoeff(&c);
                    Eigen::RowVector3d new_c = U.row(F(fid, c));
                    if (C.rows() == 0 || (C.rowwise() - new_c).rowwise().norm().minCoeff() > 0)
                    {      
                        C.conservativeResize(C.rows() + 1, 3);
                        C.row(C.rows() - 1) = new_c;
                        CI.conservativeResize(CI.rows() + 1, 1);
                        CI(CI.rows() - 1) = (int)F(fid, c);
                        set_points(C, 0);

                        vertex_added = true;
                        return true;
                    }
                }
            }
            else
            {
                if (C.rows() > 0) //if any other mouse button
                                    //we are changing whcih of the handles we are clicking
                {
                    // Move closest control point
                    Eigen::MatrixXf CP;
                    igl::project(
                        Eigen::MatrixXf(C.cast<float>()),
                        igl_v->core().view,
                        igl_v->core().proj, igl_v->core().viewport, CP);
                    Eigen::VectorXf D = (CP.rowwise() - last_mouse).rowwise().norm();
                    currI = (D.minCoeff(&currI) < 30) ? currI : -1;
                    if (currI != -1)
                    {
                        last_mouse(2) = CP(currI, 2);
                        return true;
                    }
                }
            }
            return false;
        };
        igl_v->callback_mouse_move = [&](igl::opengl::glfw::Viewer&, int, int)->bool
        {
            if (currI != -1)
            {
                Eigen::RowVector3f drag_mouse(
                    igl_v->current_mouse_x,
                    igl_v->core().viewport(3) - igl_v->current_mouse_y,
                    last_mouse(2));
                Eigen::RowVector3f drag_scene, last_scene;
                igl::unproject(
                    drag_mouse,
                    igl_v->core().view,
                    igl_v->core().proj,
                    igl_v->core().viewport,
                    drag_scene);
                igl::unproject(
                    last_mouse,
                    igl_v->core().view,
                    igl_v->core().proj,
                    igl_v->core().viewport,
                    last_scene);
                RowVector3d diff = (drag_scene - last_scene).cast<double>();
                MatrixXd Z = MatrixXd::Zero(C.rows(), C.cols());
                Z.row(currI) = diff;
                C.row(currI) += diff;
                last_mouse = drag_mouse;
                set_points(C, 0);
                return true;
            }
            return false;
        };
        igl_v->callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool
        {
            currI = -1;
            return false;
        };
	}

    VectorXi get_selected_indices()
    {
        return CI;
    }

    MatrixXd get_handle_positions()
    {
        return C;
    }

    /*
    Queries new handle information. Once this is called, the vertex_added flag is set to false
    Outputs:
        C - n x 3 vertex handle positions
        CI - n x 1 vertex handle indices
        vertex_added - (bool) whether or not a new vertex was added this timestep
    */
    bool query_new_handles(MatrixXd& C, VectorXi& CI)
    {
        C = this->C;
        CI = this->CI;

        bool new_vert = vertex_added;
        if (new_vert)
            vertex_added = false;
        return new_vert;
    }

    /*
    Before querying new handles, if a new one was added, project it to the mesh V, F
    */
    bool query_new_handles_on_mesh(MatrixXd& C, VectorXi& CI, const MatrixXd& V, const MatrixXi& F)
    {
        C = this->C;
        CI = this->CI;

        bool new_vert = vertex_added;
        if (new_vert)
        {
            VectorXd sqrD;
            MatrixXd U;
            VectorXi VI = igl::colon<int>(0, V.rows() - 1);
            igl::point_mesh_squared_distance(C, V, VI, sqrD, CI, U);
            C = U;
            this->C = C;
            this->CI = CI;
            vertex_added = false;
        }
        return new_vert;
    }

    /*
    Call this to label all the handle's indices as old
    */
    void stalify_new_vertices()
    {
        vertex_added = false;
    }
};