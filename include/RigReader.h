#pragma once
#include "rig.h"

/*
opens the .json file in rig file path and figures out what rig to use. Returns a pointer to said rig, on the heap!
*/
void rig_reader(std::string rig_file_path, Rig* rig);

/*
Same as above but also returns new cache directories associated with the rig, if required. Simply uses the default, assumes the cache is located at /rig_dir/cache/
*/
void rig_reader(std::string rig_file_path, Rig* rig, std::string& B_file_path, std::string& S_file_path, std::string& clusters_file_path);