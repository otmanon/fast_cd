#pragma once
#include <Eigen/Core>

Eigen::MatrixXd get_rainbow_colormap()
{
	//obtained from https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
	Eigen::MatrixXd rainbow_cmap = Eigen::MatrixXd(12, 3);
	//rainbow_cmap << 166, 31, 178, 51, 251, 227, 253, 255, 202, 106, 255, 177,
	//				206, 120, 223, 160, 154, 26, 191, 127, 178, 61, 255, 89,
	//				227, 180, 138, 44, 153, 28, 111, 0, 214, 154, 153, 40;

	rainbow_cmap << 166, 206, 227,
		31, 120, 180,
		178, 223, 138,
		51, 160, 44,
		251, 154, 153,
		227, 26, 28,
		253, 191, 111,
		255, 127, 0,
		202, 178, 214,
		106, 61, 154,
		255, 255, 153,
		177, 89, 40;
	rainbow_cmap /= 255;
	return rainbow_cmap;

}