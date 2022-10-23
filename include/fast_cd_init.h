#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <json.hpp>
#include <Eigen/Core>
#include <igl/list_to_matrix.h>
#include <igl/PI.h>
using namespace nlohmann;
struct FastCDInit {


    std::string json_filepath;
    json j;
    std::string mesh; //.obj mesh mpath
    std::string rig; //.json rig path for rig file
    std::string anim; //.json or .dmat anim path
    std::string rig_controller;
    
    std::string results; //results directory where we place outputs (if anya re specified)
    std::string display_mesh; //dsiplay texture .obj of the mesh with a higher resolution and UV coordinates
    std::string texture;     //texture png path 
    std::string weights;        //path to weights.DMAT file.
    std::string weights_type; //wether bbw for bounded biharmonic weights, or bh for bone heat (these should be in the rig.json file given to us from blender.)
    
    double ym, pr; //youngs modulus and poisson ratio
    double dt;    //timestep
    double beta;
    int num_substeps;
    int num_modes; int num_clusters; 
    double momentum_leaking_dt;
    int num_clustering_features;
    int substeps;
    bool do_clustering, do_reduction, do_inertia, do_cd;
    double gamma; // when 0, then no rig momentum, when 1, full rig momentum.

    double metric_alpha;
    std::string metric_type;

    std::string modes_dir;   //directory where to search for the modes. Must have a B.DMAT  for the eigenvectors and an L.DMAT file for the eigenvalues, otherwize, will be recomputed
    std::string cluster_dir; // directory for the cluster cache.

    bool screenshot; //whether or not to take screenshots. If so saves them in <results_dir>/screenshots/
    bool record_rig;
    bool record_mesh;
    bool record_modal_activations;
    bool record_energy;
    bool  record_gradient;
    bool  record_gradient_norm;
    bool record_timings;
    bool use_complementary_recordings;

    bool update_basis_naive;
    bool update_basis_fast;

    Eigen::RowVector3d color; //color to draw the mesh in, if the texture isn't found or isn't specified
    
    double rig_thickness;
    Eigen::RowVector3d rig_color;
    bool vis_rig; //if this is true we do not do complementary dynamics and ONLY visualize the rig
    bool vis_clusters;
    bool vis_texture;

    Eigen::RowVector3d eye_pos;
    Eigen::RowVector3d center;
    double zoom;
    bool face_based;

    bool init_z;
    bool randomize_init_z;
    std::string init_z_dir;
    int init_z_step;
    int init_p_step;
    bool animate_rig;

    bool run_solver_to_convergence;
    double convergence_threshold;
    double max_iters;
    /*
    sim /vis_clusters / vos_modes/ vis_weights
    */
	void init(int argc, char* argv[] )//
	{
        json_filepath = argc > 1 ? argv[1] : "../data/cluster_visualisation/charizard_null_rig.json"; ////"../data/prism_metric_tests/twisting_elastic_metric.json";
      // if (init_type == "sim")
      //     init_sim_from_json(json_filepath);
      // else if (init_type == "vis_clusters")
      //     init_cluster_vis_from_json(json_filepath);
      // else if (init_type == "vis_modes")
      //     init_mode_vis_from_json(json_filepath);
      // else if (init_type == "vis_weights")
      //     init_weight_vis_from_json(json_filepath);
      // if (init_type == "interactive")
      //     init_interactive_sim_from_json(json_filepath);
        //standard init sim... everything is not required.
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, false);
        read_json_entry_filepath(j, "rig", rig, false);
        read_json_entry_filepath(j, "anim", anim, false);

        read_json_entry_filepath(j, "display_mesh", display_mesh, false);
        read_json_entry_filepath(j, "texture", texture, false);

        read_json_entry_filepath(j, "weights", weights, false);
        weights_type = j.value("weights_type", "bbw");

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder

        read_json_entry_filepath(j, "modes_dir", modes_dir, false, fs::path(rig).parent_path().string() + "/cache/modes/default/");
        read_json_entry_filepath(j, "clusters_dir", cluster_dir, false, fs::path(rig).parent_path().string() + "/cache/clusters/default/");

        do_reduction = j.value("do_reduction", true);
        do_clustering = j.value("do_clustering", true);
        do_cd = j.value("do_cd", true);
        do_inertia = j.value("do_inertia", true);
        num_modes = j.value("num_modes", 100);
        num_clusters = j.value("num_clusters", 100);
        gamma = j.value("gamma", 1.0);


        beta = j.value("beta", 1.0);
        init_z = j.value("init_z", false);
        init_z_dir = j.value("init_z_dir", "");
        init_z_step = j.value("init_z_step", 0);
        init_p_step = j.value("init_p_step", 0);
        animate_rig = j.value("animate_rig", true);

        randomize_init_z = j.value("randomize_init_z", false);
        randomize_init_z = (randomize_init_z && !init_z); //only allow randomize if init_z unspecified

        num_clustering_features = j.value("num_clustering_features", 10);
        ym = j.value("ym", 0.1);
        pr = j.value("pr", 0.0);
        dt = j.value("dt", 1.0 / 60.0);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        metric_type = j.value("metric_type", "momentum");
         metric_alpha = j.value("metric_alpha", 1.0);

        update_basis_naive = j.value("update_basis_naive", false);
        update_basis_fast = j.value("update_basis_fast", false);

        num_substeps = j.value("num_substeps", 1);

        run_solver_to_convergence = j.value("run_solver_to_convergence", true);
        convergence_threshold = j.value("convergence_threshold", 1e-6);
        max_iters = j.value("max_iters", 100);

        screenshot = j.value("screenshot", false);
        record_rig = j.value("record_rig", false);
        record_mesh = j.value("record_mesh", false);
        record_modal_activations = j.value("record_modal_activations", false);
        record_energy = j.value("record_energy", false);
        record_gradient_norm = j.value("record_gradient_norm", false);
        record_gradient = j.value("record_gradient", false);
        record_timings = j.value("record_timings", false);

        use_complementary_recordings = j.value("use_complementary_recordings", false);
        vis_texture = j.value("vis_texture", false);

        vis_rig = j.value("vis_rig", false);
        vis_clusters = j.value("vis_clusters", false);
        rig_thickness = j.value("rig_thickness", 0.1);
        std::vector<double> baby_blue_list = { 0.537, 0.81176,  0.9411 };
        std::vector<double> dull_yellow_list = { 0.925, 0.890, 0.631 };
        std::vector<double> color_list = j.value("color", baby_blue_list);
        igl::list_to_matrix(color_list, color);

        color_list = j.value("rig_color", dull_yellow_list);
        igl::list_to_matrix(color_list, rig_color);

        //zoom
        init_cam_parameters(j);

        if (!fs::exists(fs::path(results)))
        {
            fs::create_directories(fs::path(results));
        }
        std::ofstream o(results + "init.json");
        o << std::setw(4) << j << std::endl;
	}

    void read_json_entry_filepath(json& j, std::string json_key, std::string& filepath, bool required, std::string default_path = "", bool confirm_exists=false)
    {
        namespace fs = std::filesystem;
        if (required)
        {
            if (j.count(json_key) == 0)
            {
                std::string key = json_key;
                printf("%s , did not specify required entry in init.json\n", key.c_str());
                exit(0);
            }
            filepath = j[json_key];
        }
        if (!required)
        {
            if (j.count(json_key) == 0)
            {
                filepath = default_path;
            }
            else
            {
                filepath = j[json_key];
            }
        }        
        if (!fs::exists(fs::path(filepath)) && confirm_exists)
        {
            printf("%s , could not find required file specified in init.json \n", filepath.c_str());
            exit(0);
        }
    }

    void init_cam_parameters(json j)
    {
        zoom = j.value("zoom", 1.0);
        std::vector<double> eye_pos0 = { 0, 0, 5 };
        std::vector<double> eye_pos_list = j.value("eye_pos", eye_pos0);
        igl::list_to_matrix(eye_pos_list, eye_pos);

        std::vector<double> center0 = { 0, 0, 0 };
        std::vector<double> center_list = j.value("center", center0);
        igl::list_to_matrix(center_list, center);

        face_based = j.value("face_based", true);
    }
};


struct InitSim : public FastCDInit {

    //uses default init and overrides necessary
    void init(int argc, char* argv[] )
    {
        FastCDInit::init(argc, argv );
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
        read_json_entry_filepath(j, "rig", rig, true, "", true);
        read_json_entry_filepath(j, "anim", anim, true, "", true);
    }
};

struct InitInteractiveSim : public FastCDInit {

    //uses default init and overrides necessary
    void init(int argc, char* argv[])
    {
        FastCDInit::init(argc, argv);
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
        read_json_entry_filepath(j, "rig", rig, true, "", true);
    }
};

struct InitBuildRig : public FastCDInit 
{
    std::string output_rig;


    int desired_num_bones;
    void init(int argc, char* argv[])
    {
        FastCDInit::init(argc, argv);
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
 
        read_json_entry_filepath(j, "output_rig", output_rig, true, "", false);

        desired_num_bones = j.value("desired_num_bones", 10);
    }
};


struct InitSimRigSubspace : public FastCDInit
{
    std::string output_rig;
    std::string dummy_rig;

    int desired_num_bones;
    void init(int argc, char* argv[])
    {
        FastCDInit::init(argc, argv);
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);

        read_json_entry_filepath(j, "output_rig", output_rig, true, "", false);
        read_json_entry_filepath(j, "dummy_rig", dummy_rig, true, "", true);

        desired_num_bones = j.value("desired_num_bones", 10);

        record_mesh = true;
        record_modal_activations = true;
    }
};


struct InitSimScriptedRig : public FastCDInit
{
    std::string scripted_rig_controller;

    double rad_per_frame;
    Eigen::Vector3d rot_axis;
    int max_step;
    void init(int argc, char* argv[])
    {
        
        FastCDInit::init(argc, argv);
        
        scripted_rig_controller = j.value("scripted_rig_controller", "rotater");

        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
        read_json_entry_filepath(j, "rig", rig, true, "", false);
        max_step = j.value("max_step", 360);

        ///For rotater
        double deg_per_frame = j.value("deg_per_frame", 1.0);
        rad_per_frame = deg_per_frame * igl::PI / 180.0;
        std::vector<double> axis_default = { 0, 0, 1 };
        std::vector<double> axis_list = j.value("rot_axis", axis_default);
        igl::list_to_matrix(axis_list, rot_axis);
    }
};

struct InitSimRotaterController : public  InitSimScriptedRig
{
    std::string scripted_rig_controller;

    double rad_per_frame;
    Eigen::Vector3d rot_axis;
    void init(int argc, char* argv[])
    {
        InitSimScriptedRig::init(argc, argv);


    }
};

struct InitModeVis : public FastCDInit
{
    std::string colormap;
    int num_ambient_occlusion_rays;

    bool do_ambient_occlusion;
    int period;
    double mode_scale;

    int num_cmap_steps;
    double min_cmap;
    double max_cmap;
    std::string mode_type;
    void init(int argc, char* argv[])
    {
        FastCDInit::init(argc, argv);
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
        read_json_entry_filepath(j, "rig", rig, true, "", false);

         period = j.value("period", 100);
        mode_scale = j.value("mode_scale", 1.0);

        colormap = j.value("colormap", "parula");
        num_ambient_occlusion_rays = j.value("num_ambient_occlusion_rays", 100);
        do_ambient_occlusion = j.value("do_ambient_occlusion", false);

        mode_type = j.value("mode_type", "complementary");

        num_cmap_steps = j.value("num_cmap_steps", 21);
        min_cmap = j.value("min_cmap",0.0);
        max_cmap = j.value("max_cmap", 1.0);
    }
};



struct InitClustersVis : public FastCDInit
{
  
    std::string clusters_type;
    void init(int argc, char* argv[])
    {
        FastCDInit::init(argc, argv);
        read_json_entry_filepath(j, "mesh", mesh, true, "", true);
        read_json_entry_filepath(j, "rig", rig, true, "", false);

        clusters_type = j.value("mode_type", "complementary");
    }
};