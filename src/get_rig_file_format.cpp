#include "get_rig_file_format.h"
#include <json.hpp>
#include <filesystem>
#include <fstream>

using namespace nlohmann;

std::string get_rig_file_format(std::string rig_file)
{
	
		//open surface_file as json
		namespace fs = std::filesystem;

		json j;
		std::ifstream i(rig_file);
		i >> j;

		if (j.value("format", "volume") == "volume")
		{
			return "volume";
		}
		else
		{
			return "surface";
		}

	


}