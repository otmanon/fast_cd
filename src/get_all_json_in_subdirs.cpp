#pragma once
#include "get_all_json_in_subdirs.h"

#include <filesystem>

//goes through all directories inside query_directory, Inside each subdirectory, makes sure there is a .json file. If there is, then adds that directory
//to valid directories, and adds that json filename (without the extension) to json_name.
void get_all_json_in_subdirs(const std::string& query_directory, std::vector<std::string>& valid_directories, std::vector<std::string>& json_names)
{
namespace fs = std::filesystem;


    json_names.clear();
    valid_directories.clear();
    if (!fs::exists(fs::path(query_directory)))
    {
        fs::create_directories(fs::path(query_directory));
    }
    for (const auto& entry : fs::directory_iterator(query_directory))
    {        
        if (fs::is_directory(entry))
        {
            for (const auto& sub_entry : fs::directory_iterator(entry.path()))
            {
                const std::string ext = sub_entry.path().extension().string();
                if (ext == ".json")
                {
                    json_names.push_back(sub_entry.path().stem().string());
                    valid_directories.push_back(sub_entry.path().string());
                }
            }

        }
    }
}

void get_all_json_in_dir(const std::string& query_directory, std::vector<std::string>& valid_directories, std::vector<std::string>& json_names)
{
 namespace fs = std::filesystem;


    json_names.clear();
    valid_directories.clear();
    if (!fs::exists(fs::path(query_directory)))
    {
        fs::create_directories(fs::path(query_directory));
    }
    for (const auto& entry : fs::directory_iterator(query_directory))
    {
        const std::string ext = entry.path().extension().string();
        if (ext == ".json")
        {
            json_names.push_back(entry.path().stem().string());
            valid_directories.push_back(entry.path().string());
        }
            
    }
}
