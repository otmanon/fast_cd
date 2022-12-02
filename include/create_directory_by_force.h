#pragma once
#include <filesystem>
using namespace std;

/*
* Creates a directory in filesystem, including all parent directories required from it that dont't exist.
*/
bool create_directory_by_force(string& directory)
{
	bool success = true;
	namespace fs = std::filesystem;
	if (!fs::exists(fs::path(directory)))
	{
		success = success && fs::create_directories(fs::path(directory));
	}
	return success;
}