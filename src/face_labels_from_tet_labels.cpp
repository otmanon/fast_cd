#include "face_labels_from_tet_labels.h"
Eigen::VectorXi face_labels_from_tet_labels(Eigen::VectorXi& labels,  Eigen::VectorXi& FiT)
{
	Eigen::VectorXi fC(FiT.rows());
	for (int i = 0; i < FiT.rows(); i++)
	{
		fC(i) = labels(FiT(i));
	}
	return fC;
}