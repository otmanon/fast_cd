#include "format_transforms.h"
#include <igl/PI.h>
#include <vector>

void format_transforms(const std::vector<Eigen::Affine3d>& T, std::string format, std::vector<Eigen::Affine3d>& A)
{
	A.resize(T.size());
	if (format == "blender")
	{
		for (int i = 0; i < T.size(); i++)
		{
			A[i] = T[i] * Eigen::AngleAxisd(-igl::PI * 0.5, Eigen::Vector3d::UnitZ());  // get skeleton inverse rest pos
			A[i] = Eigen::AngleAxisd(igl::PI * 0.5, Eigen::Vector3d::UnitX()); //bone resting on different axis
		}
	}
}