#include "save_rig_recording.h"
#include <igl/writeDMAT.h>

void save_rig_recording(std::string filename, Eigen::MatrixXd P)
{
	igl::writeDMAT(filename, P);
}