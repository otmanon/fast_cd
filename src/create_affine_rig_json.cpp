#include "create_affine_rig_json.h"
#include <fstream>
#include <stdio.h>
#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <ostream>
#include "json.hpp"
#include <iomanip>
void create_affine_rig_json(std::string rig_path)
{
    //write current state to json
    std::ofstream o(rig_path);
    nlohmann::json j;

    j["rig_type"] = "affine_rig";

    o << std::setw(4) << j << std::endl;
}