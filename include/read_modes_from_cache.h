#pragma once

#include "read_modes_from_cache.h"

#include <Eigen/Core>
#include <string>
#include <filesystem>
#include <json.hpp>
#include <igl/readDMAT.h>
#include <fstream>
#include <iostream>
using namespace nlohmann;
using namespace Eigen;

/*
Reads displacement modes from cache directory modes_cache_dir. Loads them up into B and L.
returns false if B or L can't be found, true if they can
*/
bool read_displacement_modes_from_cache(const std::string& modes_cache_dir,
 MatrixXd& B,  VectorXd& L)
{
    namespace fs = std::filesystem;
    B.resize(0, 0);
    L.resize(0);
    std::string dir = modes_cache_dir;
    bool found = false;
  
    
    printf("looking for cached displacement modes in %s \n", modes_cache_dir.c_str());


   found = igl::readDMAT(dir + "B.DMAT", B) && igl::readDMAT(dir + "L.DMAT", L);
       
   if (found)
   {
       assert(B.cols() == L.rows() && "Number of cached eigenvectors must be equal to number of cached eigenvalues");
       return true;
   }
   else
       return false;
}


/*
Compares two json files and specifically makes sure all parameters 
having to do with computing modes match If they don't, return false;
*/

bool check_mode_cache_safety(json& j, json& j_check)
{
    bool safe = true;
    safe = safe && j["mode_type"] == j_check["mode_type"];
    safe = safe && j["num_modes"] == j_check["num_modes"];   
}
/*
Reads displacement modes from cache directory modes_cache_dir. Loads them up into B and L.
returns false if B or L can't be found, true if they can.

In addition does a safety check to see if the modes are safe to be read. 
Checks cache_params.json in the directory modes_cache_dir and matches it to json j. If all mode-relevant quantities
match, then returns true. Otherwise returns false. If can't find the cache_params.json, returns false
*/
bool read_displacement_modes_from_cache(const std::string& modes_cache_dir, const json& j,
    MatrixXd& B, VectorXd& L)
{
    namespace fs = std::filesystem;
    B.resize(0, 0);
    L.resize(0);
    std::string dir = modes_cache_dir;
    bool found = false;


    printf("looking for cached displacement modes in %s \n", modes_cache_dir.c_str());

    std::string cache_json_file = modes_cache_dir + "cache_params.json";
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(cache_json_file)))
    {
        printf("%s , cache_params.json file not found in %s \n", (modes_cache_dir + "cache_params.json").c_str());
        return false;
    }
    else
    {
        std::ifstream i(cache_json_file);
        json cache_json;
        try
        {
            i >> cache_json;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "could not read cache_parms.json. parse error at byte " << ex.byte << std::endl;
        }
        bool safe = check_mode_cache_safety(cache_json, j);
        if (!safe)
        {
            printf("cache is not safe to read from! \n");
            return false;
        }
    }


    found = igl::readDMAT(dir + "B.DMAT", B) && igl::readDMAT(dir + "L.DMAT", L);

    if (found)
    {
        assert(B.cols() == L.rows() && "Number of cached eigenvectors must be equal to number of cached eigenvalues");
        return true;
    }
    else
        return false;
}


bool read_skinning_modes_from_cache(const std::string& modes_cache_dir,
	 MatrixXd& W,  VectorXd& L)
{
    namespace fs = std::filesystem;
    W.resize(0, 0);
    L.resize(0);
    std::string dir = modes_cache_dir;
    bool found = false;


    printf("looking for cached skinning modes in %s \n", modes_cache_dir.c_str());


    found = igl::readDMAT(dir + "W.DMAT", W) && igl::readDMAT(dir + "L.DMAT", L);

    if (found)
    {
        assert(W.cols() == L.rows() && "Number of cached eigenvectors must be equal to number of cached eigenvalues");
        return true;
    }
    else
        return false;
}



bool read_skinning_modes_from_cache(const std::string& modes_cache_dir, const json& j,
    MatrixXd& W, VectorXd& L)
{
    namespace fs = std::filesystem;
    W.resize(0, 0);
    L.resize(0);
    std::string dir = modes_cache_dir;
    bool found = false;


    printf("looking for cached displacement modes in %s \n", modes_cache_dir.c_str());

    std::string cache_json_file = modes_cache_dir + "cache_params.json";
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(cache_json_file)))
    {
        printf("%s , cache_params.json file not found in %s \n", (modes_cache_dir + "cache_params.json").c_str());
        return false;
    }
    else
    {
        std::ifstream i(cache_json_file);
        json cache_json;
        try
        {
            i >> cache_json;
        }
        catch (json::parse_error& ex)
        {
            std::cerr << "could not read cache_parms.json. parse error at byte " << ex.byte << std::endl;
        }
        bool safe = check_mode_cache_safety(cache_json, j);
        if (!safe)
        {
            printf("cache is not safe to read from! \n");
            return false;
        }
    }


    found = igl::readDMAT(dir + "W.DMAT",W) && igl::readDMAT(dir + "L.DMAT", L);

    if (found)
    {
        assert(W.cols() == L.rows() && "Number of cached eigenvectors must be equal to number of cached eigenvalues");
        return true;
    }
    else
        return false;
}