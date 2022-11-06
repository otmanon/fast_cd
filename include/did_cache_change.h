#pragma once
#include <string>
#include <fstream>
#include <filesystem>
#include <json.hpp>
#include <iostream>
using namespace nlohmann;
using namespace std;
/*
For two json files, cache_json, current_json, looks at each of their physical parameters that affect the simulation
(mode type), young's modulus etc.
If all parameters match, return false,
if one parameter doesn't match, return true
*/
bool did_cache_change(string cache_json, string current_json)
{
    namespace fs = std::filesystem;
	bool match = true;

    if (!fs::exists(fs::path(cache_json)))
    {
        //if cant find cahce json, it most definitley changed
        return true;
    }
    
    json cache_j; json current_j;

    std::ifstream i(cache_json);
    std::ifstream i2(current_json);
  
    
    try
    {
        i >> cache_j;
    }
    catch (json::parse_error& ex)
    {
        std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
    }

    try
    {
        i2 >> current_j;
    }
    catch (json::parse_error& ex)
    {
        std::cerr << "init.json parse error at byte " << ex.byte << std::endl;
    }

    match = (cache_j["ym"] == current_j["ym"]) && match;
    match = (cache_j["pr"] == current_j["pr"]) && match;
    match = (cache_j["height"] == current_j["height"]) && match;
    match = (cache_j["num_modes"] == current_j["num_modes"]) && match;
    match = (cache_j["mode_type"] == current_j["mode_type"]) && match;
    match = (cache_j["num_clusters"] == current_j["num_clusters"]) && match;
    match = (cache_j["num_clustering_features"] == current_j["num_clustering_features"]) && match;
    match = (cache_j["do_inertia"] == current_j["do_inertia"]) && match;
    match = (cache_j["dt"] == current_j["dt"]) && match;
    match = (cache_j["num_substeps"] == current_j["num_substeps"]) && match;
    match = (cache_j["do_reduction"] == current_j["do_reduction"]) && match;
    match = (cache_j["mesh"] == current_j["mesh"]) && match;
    match = (cache_j["weight"] == current_j["weight"]) && match;
    
    return !match;


}