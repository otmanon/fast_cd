#include "PhysicsHook.h"
#include "igl/eigs.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "rig.h"
#include "read_write.h"
class EigenHook : public PhysicsHook
{
public:
    EigenHook() : PhysicsHook() {}

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu)
    {

    }

    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;
    }

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer& viewer)
    {
        viewer.data().set_vertices(Q);
       

    }

    void initViewer(igl::opengl::glfw::Viewer& viewer)
    {
        this->viewer = &viewer;
        viewer.data().set_mesh(origQ, F);
    }

private:
    double k;
    double dt;
    Eigen::MatrixXd origQ, Q, V;
    Eigen::MatrixXi F, T;

    int mode = 0;
    int step = 0;
    int period = 100000;
    int r = 10;
    int half_period = period / 2;
    double scale = 1e-2;
    //eigenvectors and values of our laplacian
    Eigen::MatrixXd B_igl, B_spectra;
    Eigen::VectorXd S;

    Eigen::MatrixXd N_faces;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;

    AffineRig rig;
    igl::opengl::glfw::Viewer* viewer;
};