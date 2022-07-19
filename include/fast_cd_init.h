#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <json.hpp>
#include <Eigen/Core>
#include <igl/list_to_matrix.h>
using namespace nlohmann;
struct FastCDInit {

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
    bool do_clustering, do_reduction, do_inertia;
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

    Eigen::RowVector3d color; //color to draw the mesh in, if the texture isn't found or isn't specified
    
    double rig_thickness;
    Eigen::RowVector3d rig_color;
    bool vis_rig; //if this is true we do not do complementary dynamics and ONLY visualize the rig
    bool vis_clusters;
    bool vis_texture;

    Eigen::RowVector3d eye_pos;
    Eigen::RowVector3d center;
    double zoom;

    bool init_z;
    std::string init_z_dir;
    int init_z_step;

    /*
    sim /vis_clusters / vos_modes/ vis_weights
    */
	void init(int argc, char* argv[], std::string init_type)
	{
        std::string json_filepath = argc > 1 ? argv[1] : "../data/beta_0_equivalents/fish.json"; ////"../data/prism_metric_tests/twisting_elastic_metric.json";
        if (init_type == "sim")
            init_sim_from_json(json_filepath);
        else if (init_type == "vis_clusters")
            init_cluster_vis_from_json(json_filepath);
        else if (init_type == "vis_modes")
            init_mode_vis_from_json(json_filepath);
        else if (init_type == "vis_weights")
            init_weight_vis_from_json(json_filepath);
        if (init_type == "interactive")
            init_interactive_sim_from_json(json_filepath);

        // save this json file in the results seciton so we know what all the most up to date results are for
        
	}

    void read_json_entry_filepath(json& j, std::string json_key, std::string& filepath, bool required, std::string default = "")
    {
        namespace fs = std::filesystem;
        if (required)
        {
            if (j.count(json_key) == 0)
            {
                printf("%s , did not specify required %s entry in init.json \n", json_key.c_str());
                exit(0);
            }
            filepath = j[json_key];
        }
        if (!required)
        {
            if (j.count(json_key) == 0)
            {
                filepath = default;
            }
            else
            {
                filepath = j[json_key];
            }
        }
       
        if (!fs::exists(fs::path(filepath)) && required)
        {
            printf("%s , could not find required file specified in init.json \n", filepath.c_str());
            exit(0);
        }
    }

    void init_interactive_sim_from_json(std::string json_filepath)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        json j;
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, true);
        read_json_entry_filepath(j, "rig", rig, true);

        rig_controller = j.value("rig_controller", "handle"); // "handle", "skeletonFK", "skeletonIK"

        read_json_entry_filepath(j, "anim", anim, false);
        read_json_entry_filepath(j, "display_mesh", display_mesh, false);
        read_json_entry_filepath(j, "texture", texture, false);

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder


        read_json_entry_filepath(j, "modes_dir", modes_dir, false, fs::path(rig).parent_path().string() + "/cache/modes/default/");
        read_json_entry_filepath(j, "clusters_dir", cluster_dir, false, fs::path(rig).parent_path().string() + "/cache/clusters/default/");

        do_reduction = j.value("do_reduction", true);
        do_clustering = j.value("do_clustering", true);
        do_inertia = j.value("do_inertia", true);
        num_modes = j.value("num_modes", 100);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        num_clusters = j.value("num_clusters", 100);

        num_clustering_features = j.value("num_clustering_features", 10);
        ym = j.value("ym", 10);
        pr = j.value("pr", 0.0);
        dt = j.value("dt", 1.0 / 60.0);
        beta = j.value("beta", 1.0);
        num_substeps = j.value("num_substeps", 1);

        screenshot = j.value("screenshot", false);
        record_rig = j.value("record_rig", false);
        record_mesh = j.value("record_mesh", false);
        record_modal_activations = j.value("record_modal_activations", false);
        record_energy = j.value("record_energy", false);
        record_gradient_norm = j.value("record_gradient_norm", false);
        record_gradient = j.value("record_gradient", false);
        record_timings = j.value("record_timings", false);

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

        //write current state to json

        if (!fs::exists(fs::path(results)))
        {
            fs::create_directories(fs::path(results));
        }
        std::ofstream o(results + "init.json");
        o << std::setw(4) << j << std::endl;
    }

	void init_sim_from_json(std::string json_filepath)
	{
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        json j;
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, true);
        read_json_entry_filepath(j, "rig", rig, true);
        read_json_entry_filepath(j, "anim", anim, true);

        read_json_entry_filepath(j, "display_mesh", display_mesh, false);
        read_json_entry_filepath(j, "texture", texture, false);

        read_json_entry_filepath(j, "weights", weights, false);
        weights_type = j.value("weights_type", "bbw");

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder

        read_json_entry_filepath(j, "modes_dir", modes_dir, false, fs::path(rig).parent_path().string() + "/cache/modes/default/");
        read_json_entry_filepath(j, "clusters_dir", cluster_dir, false, fs::path(rig).parent_path().string() +"/cache/clusters/default/");

        do_reduction = j.value("do_reduction", true) ;
        do_clustering = j.value("do_clustering", true);
        do_inertia = j.value("do_inertia", true);
        num_modes = j.value("num_modes", 100);
        num_clusters = j.value("num_clusters", 100);

        beta = j.value("beta", 1.0);
        init_z = j.value("init_z", true);
        init_z_dir = j.value("init_z_dir", "");
        init_z_step = j.value("init_z_step", 0);

        num_clustering_features = j.value("num_clustering_features", 10);
        ym = j.value("ym", 0.1);
        pr = j.value("pr", 0.0);
        dt = j.value("dt", 1.0 / 60.0);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        metric_type = j.value("metric_type", "momentum");
        metric_alpha = j.value("metric_alpha", 1.0);


        num_substeps = j.value("num_substeps", 1);
      

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

    void init_cluster_vis_from_json(std::string json_filepath)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        json j;
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, true);
        read_json_entry_filepath(j, "rig", rig, true);
       

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder


        read_json_entry_filepath(j, "modes_dir", modes_dir, false, fs::path(rig).parent_path().string() + "/modes/default/");
        read_json_entry_filepath(j, "clusters_dir", cluster_dir, false, fs::path(rig).parent_path().string() + "/clusters/default/");

        do_reduction = true;
        do_clustering = true;
        do_inertia = true;
        num_modes = j.value("num_modes", 100);
        num_clusters = j.value("num_clusters", 100);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        num_clustering_features = j.value("num_clustering_features", 10);
        ym = 1;
        pr = 0;
        dt =  1.0 / 60.0;


        screenshot = j.value("screenshot", false);

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
    }

    void init_mode_vis_from_json(std::string json_filepath)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        json j;
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, true);
        read_json_entry_filepath(j, "rig", rig, true);

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder

        read_json_entry_filepath(j, "modes_dir", modes_dir, false, fs::path(rig).parent_path().string() + "/modes/default/");
        read_json_entry_filepath(j, "clusters_dir", cluster_dir, false, fs::path(rig).parent_path().string() + "/clusters/default/");

        do_reduction = true;
        do_clustering = true;
        do_inertia = true;
        num_modes = j.value("num_modes", 100);
        num_clusters = j.value("num_clusters", 100);
        momentum_leaking_dt = j.value("momentum_leaking_dt", 1e-6);
        num_clustering_features = j.value("num_clustering_features", 10);
        ym = 1;
        pr = 0;
        dt = 1.0 / 60.0;


        screenshot = j.value("screenshot", false);

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
    }

    void init_weight_vis_from_json(std::string json_filepath)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(json_filepath)))
        {
            printf("%s , initialization .json file not found \n", json_filepath.c_str());
            exit(0);
        }

        std::ifstream i(json_filepath);
        json j;
        try
        {
            i >> j;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
        }
        read_json_entry_filepath(j, "mesh", mesh, true);
        read_json_entry_filepath(j, "rig", rig, true);

        read_json_entry_filepath(j, "results", results, false, "../results/default_results/"); //where should we store the results... this should be a folder

        screenshot = j.value("screenshot", false);

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
    }
};