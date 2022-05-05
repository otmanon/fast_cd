#include "get_rig_type.h"

#include <json.hpp>
#include <iostream>

#include <stdio.h>

using namespace nlohmann;
RIG_TYPE get_rig_type(std::string rig_path, std::string& rig_type)
{
    namespace fs = std::filesystem;

    std::ifstream i(rig_path);
    json j;
    i >> j;
    rig_type = j.value("rig_type", "affine_rig");

    if (rig_type == "lbs_rig")
    {
        return RIG_TYPE::LBS_RIG;
    }
    else if (rig_type == "handle_rig")
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