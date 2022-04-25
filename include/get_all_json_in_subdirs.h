#pragma once
#include <string>
#include <vector>

//goes through all directories inside query_directory, Inside each subdirectory, makes sure there is a .json file. If there is, then adds that directory
//to valid directories, and adds that json filename (without the extension) to json_name.
void get_all_json_in_subdirs(const std::string& query_directory, std::vector<std::string>& valid_directories, std::vector<std::string>& json_names);


void get_all_json_in_dir(const std::string& query_directory, std::vector<std::string>& valid_directories, std::vector<std::string>& json_names);