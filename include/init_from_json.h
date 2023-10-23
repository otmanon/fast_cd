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
using namespace Eigen;
using namespace std;
using namespace igl;
/*
Very General json file reader. All examples either use this default json filereader to map json input to 
app parameters, or extend  it with their own code*/
struct initializer
{
    string json_filepath;  // 
    string root_dir;
    json j;
    string mesh; //.obj mesh mpath
    string rig_path; //.json rig path for rig file
    string anim_path; //.json anim path
    string anim_dir; //folder contraining json anim paths
    string rig_controller;

    string vertex_shader; //vertex_shader filepath
    string fragment_shader; //fragment_shader filepath
    int max_num_bones_vs;
    
    string results; //results directory where we place outputs (if anya re specified)
    string texture_obj; //dsiplay texture .obj of the mesh with a higher resolution and UV coordinates
    string texture_png;     //texture png path 
    string weights;        //path to weights.DMAT file.
    string weights_type; //wether bbw for bounded biharmonic weights, or bh for bone heat (these should be in the rig.json file given to us from blender.)

    double ym, pr; //youngs modulus and poisson ratio
    double dt;    //timestep

    double height;
    int num_substeps;
    int num_modes; int num_clusters;
    bool split_components;
    std::string mode_type; //displacement only for now
    bool read_cache;             //whether or not to look in cached directory before reading input

    string external_force_type;
    double external_force_magnitude;

    bool recompute_modes_if_not_found;
    bool recompute_clusters_if_not_found;
    std::string precomp_cache_dir;
    std::string modes_cache_dir;
    std::string clusters_cache_dir;
    std::string cache_dir;

    bool write_cache;   //whether or not to write to cache

    string sim_constraint_type;
    string subspace_constraint_type;
    double momentum_leaking_dt;
    double momentum_leaking_force;
    int num_clustering_features;
    int substeps;
    bool do_clustering, do_reduction, do_inertia, do_cd;
    bool do_warping;

    string constraint_type; // either "cd", "distance", or "closest_point". If "distance" 
    double constraint_radius; // the distance of threshold  if constraint type is set to distance

    string modes_dir;   //directory where to search for the modes. Must have a B.DMAT  for the eigenvectors and an L.DMAT file for the eigenvalues, otherwize, will be recomputed
    string cluster_dir; // directory for the cluster cache.

    bool screenshot; //whether or not to take screenshots. If so saves them in <results_dir>/screenshots/
    bool screenshot_rig;
    bool record_rig;
    bool record_mesh;
    bool record_modal_activations;
    bool record_kinetic_energy;
    bool record_kinetic_energy_complementary;
    bool record_energy;
    bool record_gradient;
    bool record_gradient_norm;
    bool record_timings;
    bool use_complementary_recordings;

    bool update_basis_naive;
    bool update_basis_fast;

    RowVector3d color; //color to draw the mesh in, if the texture isn't found or isn't specified

    double rig_thickness;
    RowVector3d rig_color;
    bool vis_rig; //if this is true we do not do complementary dynamics and ONLY visualize the rig
    bool vis_clusters;
    bool vis_texture;
    int max_fps;

    /*
    Camera Parameters
    */
    RowVector3d eye_pos;
    RowVector3d center;
    double zoom;
    bool face_based;
    bool invert_normals;
    bool double_sided_lighting;
    bool show_lines;
    bool show_faces; 
    bool init_z;
    bool randomize_init_z;
    bool random_amplitude;

    string init_z_dir;
    int init_z_step;
    int init_p_step;
    bool animate_rig;

    bool run_solver_to_convergence;
    double convergence_threshold;
    int max_iters;

   
    bool debug; 
    string debug_dir; 
    void init(std::string json_filepath = "../data/example_json/default_init.json") ////"../data/prism_metric_tests/twisting_elastic_metric.json";)//
    {
        this->json_filepath = json_filepath;
  
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
        root_dir = j.value("root_dir", "./"); //fs::path(json_filepath).parent_path().string();
        read_json_entry_filepath(j, "mesh", mesh, false);
        read_json_entry_filepath(j, "rig_path", rig_path, false);
        read_json_entry_filepath(j, "anim_path", anim_path, false);
        anim_dir = j.value("anim_dir", rig_path + "/../anim/");
        read_json_entry_filepath(j, "texture_obj", texture_obj, false);
        read_json_entry_filepath(j, "texture_png", texture_png, false);

        read_json_entry_filepath(j, "fragment_shader", fragment_shader, false);
        read_json_entry_filepath(j, "vertex_shader", vertex_shader, false);
        max_num_bones_vs = j.value("max_num_bones_vs", 16);
        read_json_entry_filepath(j, "results", results, false, "../results/default/"); //where should we store the results... this should be a folder

        read_json_entry_filepath(j, "cache_dir", cache_dir, false, results + "/cache/");
        read_json_entry_filepath(j, "precomp_cache_dir", precomp_cache_dir, false, cache_dir + "/precomp/");
        read_json_entry_filepath(j, "modes_cache_dir", modes_cache_dir, false, cache_dir + "/modes/");
        read_json_entry_filepath(j, "clusters_cache_dir", clusters_cache_dir, false, cache_dir + "/clusters/");
        

        debug = j.value("debug", false);
        debug_dir = j.value("debug_dir", results + "/debug/");

        read_json_entry_filepath(j, "weights", weights, false);
        weights_type = j.value("weights_type", "bbw");

        mesh = root_dir + mesh;
        rig_path = root_dir + rig_path;
        anim_path = root_dir + anim_path;
        anim_dir = root_dir + anim_dir;
        texture_obj = root_dir + texture_obj;
        texture_png = root_dir + texture_png;
        results = root_dir + results;
        cache_dir = root_dir + cache_dir;
        precomp_cache_dir = root_dir + precomp_cache_dir;
        modes_cache_dir = root_dir + modes_cache_dir;
        clusters_cache_dir = root_dir + clusters_cache_dir;
        debug_dir = root_dir + debug_dir;
        weights = root_dir + weights;
        vertex_shader = root_dir + vertex_shader;
        fragment_shader = root_dir + fragment_shader;

        do_reduction = j.value("do_reduction", true);
        do_clustering = j.value("do_clustering", true);
        do_cd = j.value("do_cd", true);
        do_inertia = j.value("do_inertia", true);
        do_warping = j.value("do_warping", false);
        sim_constraint_type = j.value("sim_constraint_type", "none");
        subspace_constraint_type = j.value("subspace_constraint_type", "cd_momentum_leak"); // can be either cd, cd_momentum_leak, or none for now
        num_modes = j.value("num_modes", 16);

        num_clusters = j.value("num_clusters", 100);
        split_components = j.value("split_components", true);
        height = j.value("height", 1.0);
        constraint_type = j.value("constraint_type", "closest_point"); // either  "distance", or "closest_point". If "distance" 
        constraint_radius = j.value("constraint_radius", 1e-2); // the distance of threshold  if constraint type is set to distance
        external_force_type = j.value("external_force_type", "none");
        external_force_magnitude  = j.value("external_force_magnitude", 1.0);

        mode_type = j.value("mode_type", "displacement");
        read_cache = j.value("read_cache", true);
        recompute_modes_if_not_found = j.value(" recompute_modes_if_not_found", true);
        recompute_clusters_if_not_found = j.value(" recompute_clusters_if_not_found", true);
        write_cache = j.value("write_cache", true);
        init_z = j.value("init_z", false);
        init_z_dir = j.value("init_z_dir", "");
        init_z_step = j.value("init_z_step", 0);
        init_p_step = j.value("init_p_step", 0);
        animate_rig = j.value("animate_rig", true);

        randomize_init_z = j.value("randomize_init_z", false);
        randomize_init_z = (randomize_init_z && !init_z); //only allow randomize if init_z unspecified
        random_amplitude = j.value("random_amplitude", 1);
        num_clustering_features = j.value("num_clustering_features", 10);
        ym = j.value("ym", 1);
        pr = j.value("pr", 0.0);
        dt = j.value("dt", 1.0 / 60.0);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        momentum_leaking_force = j.value("momentum_leaking_force", 1.0);

        update_basis_naive = j.value("update_basis_naive", false);
        update_basis_fast = j.value("update_basis_fast", false);

        num_substeps = j.value("num_substeps", 1);

        run_solver_to_convergence = j.value("run_solver_to_convergence", true);
        convergence_threshold = j.value("convergence_threshold", 1e-4);
        max_iters = j.value("max_iters", 100);

        screenshot = j.value("screenshot", false);
        screenshot_rig = j.value("screenshot_rig", false);
        record_rig = j.value("record_rig", false);
        record_mesh = j.value("record_mesh", false);
        record_modal_activations = j.value("record_modal_activations", false);
   
        record_kinetic_energy = j.value("record_kinetic_energy", false);
        record_kinetic_energy_complementary = j.value("record_kinetic_energy_complementary", false);
        
        record_energy = j.value("record_energy", false);
        record_gradient_norm = j.value("record_gradient_norm", false);
        record_gradient = j.value("record_gradient", false);
        record_timings = j.value("record_timings", false);

        use_complementary_recordings = j.value("use_complementary_recordings", false);
        vis_texture = j.value("vis_texture", false);

        max_fps = j.value("max_fps", 60);
        show_lines = j.value("show_lines", false);
        show_faces = j.value("show_faces", true);
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

    void write_init_json(std::string& path)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(path).parent_path()))
        {
            fs::create_directories(fs::path(path).parent_path());
        }
        std::ofstream o(path);
        o << std::setw(4) << j << std::endl;
    }
    void read_json_entry_filepath(json& j, std::string json_key, std::string& filepath, bool required, std::string default_path = "", bool confirm_exists = false)
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
        if (!fs::exists(fs::path(root_dir + filepath)) && confirm_exists)
        {
            printf("%s , could not find required file specified in init.json \n", (root_dir + filepath).c_str());
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
        invert_normals = j.value("invert_normals", false);
        double_sided_lighting = j.value("double_sided_lighting", false);
    }
};
