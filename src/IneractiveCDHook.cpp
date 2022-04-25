#include "InteractiveCDHook.h"

#include <stdio.h>
#include <cassert>
#include <Eigen/Sparse>


#include <filesystem>


#include "normalize_geometry_unit_height.h"
#include "igl/readMSH.h"
#include <igl/boundary_facets.h>
#include <igl/colon.h>

#include <igl/unique.h>
#include <igl/writeDMAT.h>
#include <igl/readOBJ.h>
#include <igl/slice_into.h>
#include <igl/slice.h>
#include <igl/cat.h>
#include <get_all_json_in_subdirs.h>
#include <igl/png/readPNG.h>
#include <igl/repmat.h>
#include <prolongation.h>
#include <igl/repdiag.h>

#include "SkeletonRig.h"
#include "SkeletonRigFKMouseController.h"
//TODO LBS rig and affine rig are both just a a special case of handle rig
#include "lbs_rig.h"
#include "HandleRig.h"
#include "null_rig.h"


#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>

#include "ConstraintControllerSquashStretch.h"

#include "LeftRightScriptedHandleRigController.h"
#include "HandleRigController.h"
#include "igl/barycenter.h"
#include "igl/centroid.h"
//#include "write_sparse_ijv_DMAT.h"

InteractiveCDHook::InteractiveCDHook(std::string& file, igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuizmoWidget* guizmo)
{
      namespace fs = std::filesystem;
   
    this->guizmo = guizmo;
    this->viewer = viewer;

    init_app_from_json(file);

    std::cout << "creating " << as.mesh_name << " scene...." << std::endl;
   

   // load_default_file_paths();
    init_geometry(as.mesh_file_path);


    init_rig(as.rig_file_path, as.mesh_file_path);

    init_rig_controller(rig);

    //finds all rigs found in the mesh/rigs/ directory
    get_all_json_in_subdirs(fs::path(as.mesh_file_path).parent_path().string() + "/rigs/", rig_paths, rig_names);

    init_create_rigs(fs::path(as.mesh_file_path).parent_path().string(), create_rig_names, create_rig_paths);
    //TODO: put all these in their own little struct called state and compare them against new_state;

    //init_app_state();

    //gui simulation change parameters. no longer useful now that we are on a single thread
    //new cd_sim state
    
    paused = false;
    changed_rig_type = false;


    timing_range = 120;
    timings = Eigen::VectorXd::Zero(timing_range);

    sim = FastSim(V0, T, rig->S,  as.ym, as.pr, as.dt, as.r, as.l, as.pinned_mode_dir, as.pinned_clusters_dir, as.do_reduction, as.do_clustering);
    
  
   cd_sim = FastCDSim(V0, T, rig->J, as.ym, as.pr, as.dt, as.r, as.l, as.cd_mode_dir, as.cd_clusters_dir, as.do_reduction, as.do_clustering);
   

    //get all indices in J and B.
    Eigen::VectorXi ix = ext_ind;
    Eigen::VectorXi iy = ext_ind.array() + V0.rows();
    Eigen::VectorXi iz = ext_ind.array() + V0.rows() * 2;
    iext = igl::cat(1, ix, igl::cat(1, iy, iz)); // jesus why does this concatenate with 1 = rows instead of cols

    igl::slice(cd_sim.B, iext, 1, cd_B_ext);
    igl::slice(sim.B, iext, 1, pinned_B_ext);

    igl::slice(cd_sim.J, iext, 1, cd_J_ext);		//slice J the same as B




    WB = W_low_to_high * sim.B;

    init_simulation();
    init_viewer(*viewer);


    init_modal_anim_state();

    new_as = as;
}


void InteractiveCDHook::init_vis_state()
{
    v_state.coarse_vis_id = 0;
    v_state.fine_vis_id = 1;
    v_state.vis_cd = true;
    v_state.vis_mode = VIS_MODE::MATCAP;
}

void InteractiveCDHook::init_modal_anim_state()
{
    mas.mode = 0;
    step = 0;
    mas.period = 100;
    mas.half_period = mas.period / 2;
    mas.scale = 1e-3;
}

void InteractiveCDHook::init_app_from_json(std::string& file)
{
    namespace fs = std::filesystem;

    if (fs::path("../data/scene_data/" + file).extension() != ".json")
    {
        std::cout << "please ensure there exists a default init.json file in the data/ directory";
        exit(0);
    }
    std::ifstream i(file);
    json j;
    i >> j;

    as.mesh_file_path = "../data/" + j.value("mesh_file", "");
    if (!fs::exists(fs::path(as.mesh_file_path)))
    {
        std::cout << "could not find mesh file specified in init.json" << std::endl;
        exit(0);
    }
    as.mesh_name = fs::path(as.mesh_file_path).stem().string();
    as.animation_name = j.value("animation_name", as.mesh_name);
    //mesh_name = mesh_name.substr(0, mesh_name.size() - 5);


    std::string constraint_type_str = j.value("constraint_type", "complementary_dynamics");  //should be complementary dynamics or pineed
    if (constraint_type_str == "pinned") as.constraint_type = CONSTRAINT_TYPE::PINNING;
    if (constraint_type_str == "complementary_dynamics") as.constraint_type = CONSTRAINT_TYPE::COMPLEMENTARY_DYNAMICS;
    as.do_reduction = j["use_reduction"] ;
    as.do_clustering = j["use_clustering"];
    as.r = j["r"];
    as.l = j["l"];
    as.feat = j.value("clustering_modal_features", 10);
    as.ym = j["youngs_modulus"];
    as.pr = j["poisson_ratio"];

    as.dt = j["dt"]; // 1.0 / 60;

    //what kind of rig are we using? default is affine rig for now
    std::string rig_path = j.value("rig_file", "scene_data/" + as.mesh_name + "rigs/affine_rig/affine_rig.json");
    as.rig_file_path = "../data/" + rig_path;
    as.rig_controller_name = j.value("rig_controller", "HandleRigMouseController");
    as.display_file_path =  "../data/" + j.value("display_mesh", fs::path(as.mesh_file_path).parent_path().string() + "/" + as.mesh_name + ".obj");

    // hold off on doing this for now until we know how to also scale down the rig parameters too
    as.texture_png_path = "../data/" + j.value("texture_png", "");

    as.record_metrics = j.value("record_metrics", false);           //should we record simulation metrics like energy and time?
    as.record_timelength = j.value("record_timelength", 100);       //if above is true, how long is the increment we should record for? default 100
    as.results_dir = "../data/" + j.value("results_dir", "results/" + as.animation_name + "/"); //where should we store the results... this should be a folder
    as.matcap_file = j["matcap_file"];
    //animation mode: floaty play around, visualize modes
    as.animation_mode = j["animation_mode"];

    
    as.cd_mode_dir = "../data/" + j.value("cd_modes_dir", fs::path(rig_path).parent_path().string() + "/cache/modes/cd_default/"); //should get some for cd clusters and for constrained clusters
    as.cd_clusters_dir = "../data/" + j.value("cd_clusters_dir", fs::path(rig_path).parent_path().string() + "/cache/clusters/cd_default/" );

    as.pinned_mode_dir = "../data/" + j.value("pinned_modes_dir", fs::path(rig_path).parent_path().string() + "/cache/modes/pinned_default/"); //should get some for cd clusters and for constrained clusters
    as.pinned_clusters_dir = "../data/" + j.value("pinned_clusters_dir", fs::path(rig_path).parent_path().string() + "/cache/clusters/pinned_default/");
    //set this to true always for now
    as.use_inertia = true;

    new_as = as;
}


void InteractiveCDHook::init_rig(std::string& rig_file, std::string& mesh_filepath)
{
    //figure out rig type
    namespace fs = std::filesystem;
  
    if (!fs::exists(fs::path(rig_file)))
    {
        std::cout << "Could not find rig file, please provide a valid rig.json file" << std::endl;
        exit(0);
    }

    std::ifstream i(rig_file);
    json j;
    i >> j;

    std::string rig_type_str;
    as.rig_type = get_rig_type(rig_file, rig_type_str);

    bool rig_file_exists = fs::exists(fs::path(rig_file)); //assume this rig file exists otherwise we quite the program.
  
    std::cout << "Using " << rig_type_str << " rig type..." << std::endl;
    
    std::string rig_file_name = fs::path(rig_file).stem().string();
    std::string cache_dir =  fs::path(mesh_filepath).parent_path().string() + "/rigs/" + rig_file_name+ "/cache/";     

    as.rig_id = -1;

    new_as = as;

    //build our affine rig ... shouldn't do this here, TODO: make a rig builder class that initializes the rig for us. 
    //We should be building the rig in the same spot we call load_rig_file_paths
    std::string rig_file_format;
    switch (as.rig_type)
    {
    case RIG_TYPE::LBS_RIG:
        // check if the LBS_rig is in surface format.
        rig_file_format = get_rig_file_format(rig_file);
        if (rig_file_format == "surface")
        {
            rig = new SkeletonRig(rig_file, V0, T);
        }
        else
        {
            rig = new SkeletonRig(rig_file);
        }
      
        break;
    case RIG_TYPE::AFFINE_RIG:
        rig = new HandleRig(V0, T);
       // as.constraint_controller = new ConstraintController(viewer, guizmo); //null constraints
        break;
    case RIG_TYPE::NULL_RIG:
        rig = new NullRig(V0, T); 
        std::cout << "" << std::endl;
        break;
    case RIG_TYPE::HANDLE_RIG:
        rig = new HandleRig(rig_file);
      //  as.constraint_controller = new ConstraintController(viewer, guizmo);
       // exit(0);
        break;
    }

}

void InteractiveCDHook::init_rig_controller(Rig * rig)
{
    if (as.rig_controller_name == "HandleRigMouseController")
    {
        as.rig_controller = new HandleRigMouseController(rig->p0, viewer, guizmo);
    }
    else if (as.rig_controller_name == "LeftRightScriptedHandleRigController")
    {
        as.rig_controller = new LeftRightScriptedHandleRigController(rig->p0, viewer, guizmo);
    }
    else if (as.rig_controller_name == "SkeletonRigFKMouseController")
    {
        //make sure the rig is of a Skeleton type
        //only give rest parameters...
       as.rig_controller = new SkeletonRigFKMouseController(rig->p0, ((SkeletonRig *)rig)->pI, ((SkeletonRig*)rig)->lengths, viewer, guizmo);
    }
    //TODO add one more that is for IK
}
void InteractiveCDHook::init_constraint(FastCDSim& cd_sim)
{
    
    //sets default constraints based on rig types... TODO: Should be able to override this in the init.json file, 
 
}


void InteractiveCDHook::init_geometry(std::string& mesh_file)
{
    Eigen::VectorXi tritag, tettag;
    Eigen::MatrixXd V_orig;
    //read tet mesh. F is filled with junk after this
    if (!igl::readMSH(as.mesh_file_path, V_orig, F, T, tritag, tettag))
    {
        std::cout << "Could not find mesh " + as.mesh_file_path << std::endl;
    }
    as.l = as.do_clustering ? as.l : T.rows();           //if we don't do clustering, then set num clusters to num tets
    Eigen::VectorXi _n; //throwaway var
    igl::boundary_facets(T, F, FiT, _n);
    //get list of exterior vertex indices
    igl::unique(F, ext_ind);


    Eigen::MatrixXd N;  Eigen::MatrixXi FN; Eigen::MatrixXd V_h;
    has_display_mesh = igl::readOBJ(as.display_file_path, V_high_res, UV_high_res, N, F_high_res, FUV_high_res, FN);
    if (has_display_mesh)
    {
        Eigen::SparseMatrix<double> W;
        //V_high_res *= 0.99;l
        prolongation(V_high_res, V_orig, T, W);
        W_low_to_high =igl::repdiag(W, 3);
    }
    else
    {
        V_high_res = V_orig;        //simulation mesh is display mesh!
        F_high_res = F;
        W_low_to_high.resize(3*V_orig.rows(), 3*V_orig.rows());
        W_low_to_high.setIdentity();
    }

    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    V0 = normalize_geometry_unit_height(V_orig);  //scale uniformly to unit height
    V_high_res = normalize_geometry_unit_height(V_high_res);
   
    Eigen::RowVector3d t;
    igl::centroid(V0, F, t);
    V0.array().rowwise() -=  t.array(); //center geometry 

    V_high_res.array().rowwise() -= t.array();


    V_high_res0 = V_high_res;
    V = V0;

    //build F_ext, convert F to index V_ext, which contains only surface vertices
    Eigen::VectorXi ind_v;      //contains the opposite info of ext_ind, given a full volumetric vertex list TV, finds surface Verts FV
    Eigen::VectorXi a;
    igl::colon(0, ext_ind.rows() - 1, a);
    ind_v.resize(V.rows(), 1);
    ind_v.setConstant(-1);
    igl::slice_into(a, ext_ind, ind_v);

    F_ext.resizeLike(F);
    for (int i = 0; i < F.rows(); i++)
    {
        F_ext(i, 0) = ind_v(F(i, 0));
        F_ext(i, 1) = ind_v(F(i, 1));
        F_ext(i, 2) = ind_v(F(i, 2));
    }

}


//physics related
void InteractiveCDHook::init_simulation()
{
    //TODO: reset render
    V = V0;         //V tracks the deforming mesh
    igl::slice(V, ext_ind, 1, V_ext);
    V0_ext = V_ext;
    uc_curr = Eigen::VectorXd::Zero(V.rows() * V.cols()); uc_prev = uc_curr;
    u_curr = Eigen::VectorXd::Zero(V.rows() * V.cols());
    u_prev = Eigen::VectorXd::Zero(V.rows() * V.cols());
    z_prev = Eigen::VectorXd::Zero(as.r); z_curr = z_prev;

    //reset rig
    rig->reset();
    as.rig_controller->reset();
  //  as.constraint_controller->reset();
    p_curr = as.rig_controller->p_rel; p_prev=  as.rig_controller->p_rel;

    //modal animation reset
    init_modal_anim_state();
}

void InteractiveCDHook::init_create_rigs(std::string mesh_dir, std::vector<std::string>& create_rig_names, std::vector<std::string>& create_rig_paths)
{
    std::string two_handle_name = "two_handle_rig";
    std::string two_handle_filepath = mesh_dir + "/rigs/" + two_handle_name + "/" + two_handle_name + ".json";

    //create_two_handle_rig(two_handle_filepath, V);
    create_rig_paths.push_back(two_handle_filepath);
    create_rig_names.push_back(two_handle_name);

}
