#pragma once
//#include "PhysicsHook.h"
#include "igl/opengl/gl.h"
#include "lbs_rig.h"
#include "igl/min_quad_with_fixed.h"
#include "fast_complementary_dynamics_sim.h"
#include "fast_sim.h"
#include "RigController.h"
//#include "ConstraintController.h" //TODO make this its own rig tbh that way we can use model reduction and have all modes be orthogonal to the constraints
enum ANIMATION_MODE{ INTERACTIVE_ANIMATION, EIGENMODES_ANIMATION};
enum CONSTRAINT_TYPE{COMPLEMENTARY_DYNAMICS, PINNING};
enum VIS_MODE { CLUSTERS, MATCAP, TEXTURES };
class InteractiveCDHook
{
public:

    // Initializes with a file. Either a .init json file, or a mesh name
    InteractiveCDHook(std::string file, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo);

    void init_app_from_json(std::string& file);

    void init_modal_anim_state();

    void init_rig(std::string& rig_file, std::string& mesh_name);

    void init_rig_controller(Rig* rig);
    void pick_rig_controller(Rig* rig);
    void init_constraint(FastCDSim& cd_sim);
    
    void init_create_rigs(std::string mesh_dir,std::vector<std::string>& create_rig_paths, std::vector<std::string>& create_rig_names);

    void init_geometry(std::string& mesh_file);

    void init_vis_state();

    /*
    Initializes simulation parameters
    */
    virtual void init_simulation();

    /*
    Performs one typestep. Split up into two options, one if we want to look at the modes, the other if we wanna do physics. TODO: split these into two further groups
    */
    virtual bool simulateOneStep();
 
    void render_full(igl::opengl::glfw::Viewer& v, Eigen::MatrixXd& X, int cid);

    void render_reduced_cpu_proj(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::VectorXf& p, Eigen::MatrixXd& B, Eigen::SparseMatrix<double>& J, int cid);

    void render_reduced_gpu_proj(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::VectorXf & p, int n, int cid);

    void render_reduced_cpu_proj_pin(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::MatrixXd& B, Eigen::MatrixXd& X, int cid);

    void render_reduced_gpu_proj_pin(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, int n, int cid);

    void full_sim_step_pinning_control();

    void reduced_sim_step_pinning_control();

    void full_sim_step_cd_control();

    void reduced_sim_step_cd_control();

    void sim_step_modal_animation();

    //render/visualization/interaction
    virtual bool render(igl::opengl::glfw::Viewer& viewer);

    void render_cd(igl::opengl::glfw::Viewer& viewer);

    void render_pinning(igl::opengl::glfw::Viewer& viewer);
    virtual void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);

    void init_viewer(igl::opengl::glfw::Viewer& viewer);

    void init_gizmo(igl::opengl::glfw::imgui::ImGuizmoWidget& guizmo)
    {
        this->guizmo = &guizmo;
      //  rig->init_gizmo(this->guizmo);
    }

    void set_viewer_matcap(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F,  std::string matcap_file, int cid=0, int fid=1);

    void set_viewer_clusters(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters, int cid = 0, int fid = 1);

    void set_viewer_color_textures(igl::opengl::glfw::Viewer& viewer, std::string texture_filepath, Eigen::MatrixXd& V_coarse, Eigen::MatrixXi& F_coarse, Eigen::MatrixXd& V_fine, Eigen::MatrixXi& F_fine,
        Eigen::MatrixXd& UV_fine, Eigen::MatrixXi& FUV_fine, int cid, int fid);

    void set_viewer_defo_textures(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& X, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W, int cid=0, int fid=1);

    void poll_sim_changes();

    //should have a rigReader class that takes this in as an option
    void change_rig_type();
    
    void change_constraint(FastCDSim& cd_sim, CONSTRAINT_TYPE constraint_type);

    bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

    bool mouse_move(igl::opengl::glfw::Viewer& viewer, int x, int y);

    bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

    bool key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier);

    void change_animation_mode();

    //system
    void save_params();

    void save_results();
    RIG_TYPE get_rig_type(std::string rig_path, std::string& rig_type);

    void save_params(std::string custom_filename);
    //reduced_solver functions


public:


    std::vector<std::string> rig_paths, rig_names, create_rig_names, create_rig_paths;
    struct app_state {
        //parameters that will frequently be changed.
        //young's modulus and poisson ratio
        float ym, pr;
        //timestep size
        double dt;
        //number of modes
        int r;
        //number of clusters
        int l;

        //number of features used for nodal features
        int feat; 

        //turn on inertia
        bool use_inertia;
        //whether to use model reduction
        bool do_reduction;
        //whether to use cluserting
        bool do_clustering;

        //whether we should record our energy/time infor
        bool record_metrics;
        int record_timelength;       
        std::string results_dir;

        //what kind of rig are we using? default is affine rig for now
        RIG_TYPE rig_type; // we dont need this
        std::string rig_controller_name;

        CONSTRAINT_TYPE constraint_type;
        ANIMATION_MODE animation_mode;

        //simulation relevant strings, filedirs, names etc
        std::string mesh_name, animation_name, rig_anim_dir;
        std::string cd_mode_dir, cd_clusters_dir, labels_file_path, rig_file_path,
            pinned_mode_dir,
            pinned_clusters_dir;
        std::string mesh_file_path, matcap_file;
        std::string display_file_path, texture_png_path;
        int rig_id;
        int proj_gpu;
        
        float k; //TODO remove this and replace with young's modulus/poisson_ratio

        RigController* rig_controller;
        
    } as;
    app_state new_as; //new app state at each frame used to check if the app state has changed

    
    Eigen::VectorXd  uc_curr, uc_prev;  //flattened displacement, rest position, current displacement and prev displacement for inertia 
    Eigen::VectorXd  u_curr, u_prev;
    Eigen::VectorXd z_next, z_prev, z_curr;                  //reduced coefficients
    Eigen::VectorXd p_next, p_prev, p_curr;
    
    int max_iters;                                    //max number of local-global iterations                     
    double tol;                                            //convergence tolerance of local-global iterations

    // rendering parameters. none of these are directly used in the simulation, just in the display
    //mesh parameters


    struct vis_state
    {
        int coarse_vis_id;      //two data entreis, one for coarse mesh, one for fines
        int fine_vis_id;        //two data entries, one for coarse, one for fine
        bool vis_cd;
        bool show_cage; // if vis_mode==Textures, whether or not to render the cage with it
        VIS_MODE vis_mode;
    } v_state, new_v_state;

    Eigen::Matrix< float, Eigen::Dynamic, 3, Eigen::RowMajor> tex;
    Eigen::MatrixXd time_energy;       //stores result data during a simulation, and saves it to disk when the simulation is over.

    Eigen::MatrixXd V, V0;      //mesh geometry (deformed) , original mesh geometry
    Eigen::MatrixXi F_ext;      // Face indices of exterior vertices V_ext. 

    Eigen::MatrixXd V_ext, V0_ext;      //traingle mesh vertices of coarse simulation mesh (does not included interior tet mesh verts) V_ext = V[ext_ind, :]... this is only for the viewer, not used anywhere in sim. 
    
    bool has_display_mesh;                  //could we find our display mesh? 
    Eigen::MatrixXd V_high_res, V_high_res0; Eigen::MatrixXi F_high_res; // triangle mesh that is embedded inside the coarse outer mesh.  Only used in viewer for high res visualization with textures.
    Eigen::MatrixXd UV_high_res; Eigen::MatrixXi FUV_high_res; // UV coordinates and face indices into these coordinates... used for texturing.;

    Eigen::SparseMatrix<double> W_low_to_high;          //matrix that takes the low dimensional coarse outer mesh values and projects them to high dimensinoal mesh values. V_high_res = W * V

    Eigen::MatrixXd cd_B_ext, pinned_B_ext;
    Eigen::SparseMatrix<double> cd_J_ext;

    Eigen::SparseMatrix<double> cd_WJ;
    Eigen::MatrixXd cd_WB, pinned_WB; // may need to precompute this further...
    Eigen::VectorXi ext_ind, iext;    //list of exterior vertex indices in our list. iext is the same as ext_ind, but indexes all entrices in a flattened V matrix.

    Eigen::MatrixXi F, T;       //faces and tets (these are shared with the simulation... shouldn't change ever)
    Eigen::VectorXi FiT;        //Faces in Tets... which tet index does each face in F belong to?



    //rig should be able to support affine rig and lbs rig
    Rig* rig;

    int step; //what simulation time step are we on?

    //animation parameters that deform mesh to show modes, modal animation state. 
    struct mode_anim_state {
        int mode;
        int period;
        int half_period;
        double scale;
    } mas;

    //keeps track of frame timing
    int timing_range;
    Eigen::VectorXd timings;


    bool changed_rig_type;
    bool refresh;
    bool paused;

    RIG_TYPE new_rig_type;

//   bool new_do_reduction;
    //viewer
    igl::opengl::glfw::Viewer* viewer;
    igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo;
    FastCDSim cd_sim;
    FastSim sim;

 
};