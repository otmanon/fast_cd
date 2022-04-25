#include "create_affine_rig_json.h"
#include <fstream>
#include <stdio.h>
#include <json.hpp>
#include <filesystem>
#include <ostream>
void create_affine_rig_json(std::string rig_path)
{
    //write current state to json
    std::ofstream o(rig_path);
    nlohmann::json j;

    j["rig_type"] = "affine_rig";

    o << std::setw(4) << j << std::endl;
}