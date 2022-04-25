#include "InteractiveCDHook.h"
#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <stdio.h>
#include <cassert>
#include <json.hpp>
#include <create_affine_rig_json.h>
#include <get_all_json_in_subdirs.h>
#include <igl/writeDMAT.h>

void InteractiveCDHook::save_params(std::string custom_init_file)
{
    /*
#ifdef WIN32
    namespace fs = std::filesystem;
#else
    namespace fs = std::experimental::filesystem;
#endif

    //write current state to json
    std::ofstream o(custom_init_file + ".json");
    json j;

    std::string mesh_file = fs::relative(fs::path(as.mesh_file_path), fs::path("../data/scene_data/")).string();
    //in case on windows:
    std::replace(mesh_file.begin(), mesh_file.end(), '\\', '/');
    j["mesh_file"] = mesh_file;
    j["matcap_file"] = as.matcap_file;
    j["use_reduction"] = as.do_reduction;
    j["use_clustering"] = as.do_clustering;
    j["r"] = as.r;
    j["l"] = as.l;
    j["dt"] = as.dt;
    j["stiffness"] = as.k;

    //what kind of rig are we using? default is affine rig for now
    j["rig_type"] = as.rig_type;
    j["animation_mode"] = as.animation_mode;
    std::string rig_file = fs::relative(fs::path(as.rig_file_path), fs::path("../data/scene_data/")).string();
    std::replace(rig_file.begin(), rig_file.end(), '\\', '/');
    j["rig_file_dir"] = rig_file;

    o << std::setw(4) << j << std::endl;*/
}

void InteractiveCDHook::save_results()
{
    #ifdef WIN32
        namespace fs = std::filesystem;
    #else
        namespace fs = std::experimental::filesystem;
    #endif
    if (!fs::exists(as.results_dir))
    {
        fs::create_directories(fs::path(as.results_dir));
    }

    if (as.record_metrics)
    {
        igl::writeDMAT(as.results_dir + "energy_timings.DMAT", time_energy);
    }
}


RIG_TYPE InteractiveCDHook::get_rig_type(std::string rig_path, std::string& rig_type)
{
#ifdef WIN32
    namespace fs = std::filesystem;
#else
    namespace fs = std::experimental::filesystem;
#endif
    std::ifstream i(rig_path);
    json j;
    i >> j;
    rig_type = j.value("rig_type", "affine_rig");

    if (rig_type == "lbs_rig")
    {
        return RIG_TYPE::LBS_RIG;
    }
    else if(rig_type == "handle_rig")
    {
        return RIG_TYPE::HANDLE_RIG;
    }
    else if (rig_type == "null_rig")
    {
        return RIG_TYPE::NULL_RIG;
    }
    else
    {
        return RIG_TYPE::AFFINE_RIG;
    }
 

}
