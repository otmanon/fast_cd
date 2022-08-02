#pragma once
#include <igl/colormap.h>
bool get_colormap_type(std::string s, igl::ColorMapType& type)
{

	if (s == "parula")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_PARULA;
		return true;
	}
	else if (s == "jet")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_JET;
		return true;
	}
	else if (s == "plasma")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_PLASMA;
		return true;
	}
	else if (s == "magma")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_MAGMA;
		return true;
	}
	else if (s == "turbo")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_TURBO;
		return true;
	}
	else if (s == "viridis")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_VIRIDIS;
		return true;
	}
	else if (s == "inferno")
	{
		type = igl::ColorMapType::COLOR_MAP_TYPE_INFERNO;
		return true;
	}
	else 
	{
		return false;
	}
}