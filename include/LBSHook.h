#pragma once
#include "PhysicsHook.h"
#include "lbs_rig.h"
#include "igl/min_quad_with_fixed.h"
class LBSHook : public PhysicsHook
{
public:
    LBSHook(std::string mesh_name);

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu);

    virtual void initSimulation();

    /*
    Initializes common simulation state variables like uc, ur, x, u_prev, u_curr, z_prev, z_curr, z
    */
    void LBSHook::init_simulation_state();

    virtual void updateRenderGeometry() {};

    void makeClusters();

    void fillModes(const Eigen::SparseMatrix<double>& H, const Eigen::SparseMatrix<double> M, Rig& rig, int num_modes);

    void compute_grouping_matrix(int num_modes, int num_clusters);

    virtual bool simulateOneStep();

    void local_global_solver_reduced(Eigen::VectorXd& z, Eigen::VectorXd& z_next);

    void local_step_reduced(Eigen::VectorXd& z, Eigen::MatrixXd& R);

    void global_step_reduced(Eigen::MatrixXd& R, Eigen::VectorXd& z_next);

    void local_global_solver(Eigen::MatrixXd& U, Eigen::MatrixXd& U_next);

    void local_step(Eigen::VectorXd& u, Eigen::MatrixXd& R);

    void global_step(Eigen::MatrixXd& R, Eigen::VectorXd& u);

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer& viewer);


    void initViewer(igl::opengl::glfw::Viewer& viewer);


    void load_file_paths();

    void compute_cd_equality_constraints();

    void change_stiffness();

    void build_A_matrix();

    void build_CSM();

    void build_reduced_matrices();

    void poll_sim_changes();


    bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

    bool mouse_move(igl::opengl::glfw::Viewer& viewer, int x, int y);

    bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

public:

    //simulation relevant strings, filedirs, names etc
    std::string mesh_name;
    std::string data_dir;
    std::string cache_dir;
    std::string B_file_dir, S_file_dir, labels_file_dir, rig_file_dir;


    //parameters that will frequently be changed.
    //stiffness parameter
    float k;
    //timestep size
    double dt;
    //number of modes
    int r;
    //number of clusters
    int l;
    //turn on inertia
    bool use_inertia;
    //use complementary dynamics (with affine rig) TODO:generalize this with enum across multiple rigs
    bool use_cd;
    //whether to use model reduction
    bool do_reduction;
    //whether to use cluserting
    bool do_clustering;

    //simulation variables:
    Eigen::VectorXi ui;  //constrained vertex indices of ui . Basically just for each handle [handleI*3, handleI*3 + 1, handleI*3 + 2]
    Eigen::VectorXd ub;   //constrained vertex displacements. 

    //constraints change booleans
    bool changed_constraint;            //master change. set this to true if any of the other below are true as well
    bool moved_handle;                  // we moved a handle!
    bool changed_num_handles;           //we added/removed a handle (for now only adding is supported)

    Eigen::MatrixXd U, X, Ur, Uc;               // displacement, rest position
    Eigen::VectorXd u, x, u_prev, u_curr;  //flattened displacement, rest position, current displacement and prev displacement for inertia 
    Eigen::VectorXd uc, ur;                //flattened complementary motion, flattened rig displacement motion  
  
    Eigen::VectorXd z;                  //reduced coefficients
    //reduced_solver parameters
    igl::min_quad_with_fixed_data<double> precomp;          //prefactorization for linear equality constraint solvers
    Eigen::LLT<Eigen::MatrixXd> precomp_reduced;               //prefactorization for reduced linear equality constraint solvers
    int max_iters;                                    //max number of local-global iterations                     
    int tol;                                            //convergence tolerance of local-global iterations


    // rendering parameters. none of these are directly used in the simulation, just in the display
    Eigen::MatrixXd handleV;    //handle positions
    Eigen::VectorXi handleI;    //indices into V that are constrained vertices
    int currentHandle;          //which handle index into V are we currently moving
    bool mouse_dragging;        //is mouse dragging

    //mesh parameters
    Eigen::MatrixXd V, V0;      //mesh geometry (deformed) , original mesh geometry
    Eigen::MatrixXd V_ext;      //traingle mesh vertices (does not included interior tet mesh verts)
    Eigen::MatrixXd B_ext;      //B, but with only the slices corresponding to exterier vertices taken out. 
    Eigen::VectorXi ext_ind;    //list of exterior vertex indices in our list. 

    Eigen::MatrixXi F, T;       //faces and tets (these are shared with the simulation... shouldn't change ever)
    Eigen::VectorXi FiT;        //Faces in Tets... which tet index does each face in F belong to?
    //Important matrices: ARAP Hessian, cotan-Laplacian, CSM, massmatrix
    //H - arap hessian
    //C - cotan laplacian
    //CSM - covariance scatter matrix 
    //M - mass matrix
    //A - current system matrix. 
    //A_static -> quasisteady state with no inertia
    //Q -reduced version of A
    //CSMB - reduced version of CSM*B
    //CSMr - precomputed each timestep CSM*(ur + x)
    Eigen::SparseMatrix<double> H, C, CSM, M, A, A_static;

    Eigen::MatrixXd Q, CSMB, CSMr;
    Eigen::MatrixXd BTA, BTM;

    // 3*m x 3*|T| grouping matrix used to cluster two or more tets together
    Eigen::SparseMatrix<double> G;
    Eigen::VectorXi labels;
    Eigen::VectorXd linear_term_reduced;


    Eigen::SparseMatrix<double> Aeq;    //linear equality constraint lhs term
    Eigen::VectorXd Beq;                //linear equality constrain rhs term
    Eigen::MatrixXd Qeq;                //linear equality constraint on reduced system lhs term
    Eigen::VectorXd qeq;

    //rig
    LBSRig rig;

    //animation parameters that deform mesh to show modes. 
    int mode = 0;
    int step = 0;
    int period = 100000;
    int half_period = period / 2;
    double scale = 1e-2;
    int timing_range;
    Eigen::VectorXd timings;

    //reduction parameters...eigenvectors and values of our hessian
    Eigen::MatrixXd B_full;
    Eigen::VectorXd S_full;
    Eigen::MatrixXd B;
    Eigen::VectorXd S;

    bool changed_stiffness;     //turns on in render thread GUI so that we update the stiffness in the next simulation reduced_step
    float new_stiffness;
    bool switched_clustering;      //if user clicked on checkpox that switches reduced mode on/off
    bool new_do_clustering;
    int new_num_modes;
    int new_num_clusters;
    //   bool new_do_reduction;
        //viewer
    igl::opengl::glfw::Viewer* viewer;
    Eigen::Vector3d mouse_win, mouse_world;
    Eigen::Vector3d mouse_drag_win, mouse_drag_world;
};