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
Reads clusters from cache directory modes_cache_dir. Loads them up into B and L.
returns false if B or L can't be found, true if they can
*/
bool read_clusters_from_cache(const std::string& clusters_cache_dir,
    VectorXi& l)
{
    namespace fs = std::filesystem;
    l.resize(0);
    std::string dir = clusters_cache_dir;
    bool found = false;


    printf("looking for cached clusters in %s \n", clusters_cache_dir.c_str());


    found = igl::readDMAT(dir + "l.DMAT", l);

    return found;
}
