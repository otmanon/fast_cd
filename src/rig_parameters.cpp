#include "rig_parameters.h"
#include <Eigen/Dense>
/*
Translates _stacked_ affine rig parameters in P0 all by t
*/
void translate_rig_parameters(MatrixXd& P0, Vector3d& t)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	int num_b = P0.rows() / 4;
	for (int i = 0; i < num_b; i++)
	{
		P0.row(4 * i) = t.transpose();
	}
};



/*
Transforms _stacked_ affine rig parameters in P0 all by T
*/
void transform_rig_parameters(MatrixXd& P0, MatrixXd& T)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	assert((T.rows() == 4 && T.cols() == 3) || (T.rows() == 3 && T.cols() == 4) && "expecting T to be 4 x 3 matrix or a 3x 4 matrix");
	if (T.rows() == 4 && T.cols() == 3)
	{
		T = T.transpose().eval();
	}
	int num_b = P0.rows() / 4;
	Matrix4d A;
	for (int i = 0; i < num_b; i++)
	{
		A.setIdentity();
		A.block(0, 0, 3, 4) = P0.block(4 * i, 0, 4, 3).transpose();
		Eigen::MatrixXd B = T * A;
		P0.block(4 * i, 0, 4, 3) = (B).transpose();
	}
};

void transform_rig_parameters(VectorXd& p, MatrixXd& A)
{
	assert((p.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p.rows() / 12;
	int dim = 3;
	MatrixXd P = Map<MatrixXd>(p.data(), num_b * (dim + 1), dim);
	transform_rig_parameters(P, A);
	p = Map<VectorXd>(P.data(), (dim * (dim + 1)) * num_b);
}

/*
Transforms a rig animation. Given a stacked list of flattened rig parameters, apply the transform A to them
*/
void transform_rig_parameters_anim(MatrixXd& anim_P, MatrixXd& A)
{
	assert((A.rows() == 4 && A.cols() == 3) || (A.rows() == 3 && A.cols() == 4) && "expecting T to be 4 x 3 matrix or a 3x 4 matrix");


	for (int i = 0; i < anim_P.cols(); i++)
	{
		Eigen::VectorXd p = anim_P.col(i);
		transform_rig_parameters(p, A);
		anim_P.col(i) = p;
	}
}

/*
Translates _stacked_ affine rig parameters in P0 all by t on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void rotate_rig_parameters(VectorXd& p0, Vector3d& t)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P0 = Map<MatrixXd>(p0.data(), num_b * (dim + 1), dim);
	translate_rig_parameters(P0, t);
	p0 = Map<VectorXd>(P0.data(), (dim * (dim + 1)) * num_b);
}

/*
Rotates _stacked_ affine rig parameters in P0 all by R
*/
void rotate_rig_parameters(MatrixXd& P0, Matrix3d& R)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	int num_b = P0.rows() / 4;
	for (int i = 0; i < num_b; i++)
	{
		P0.block(4 * i, 0, 4, 3) = (R * P0.block(4 * i, 0, 4, 3).transpose()).transpose();
	}
};


/*
Rotates _stacked_ affine rig parameters in P0 all by R on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void rotate_rig_parameters(VectorXd& p0, Matrix3d& R)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P0 = Map<MatrixXd>(p0.data(), num_b * (dim + 1), dim);
	rotate_rig_parameters(P0, R);
	p0 = Map<VectorXd>(P0.data(), (dim * (dim + 1)) * num_b);
}

/*
Sets all rig parameters to the same affine transform T
*/
void set_rig_parameters(MatrixXd& P0, MatrixXd& T)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	assert((T.rows() == 4 && T.cols() == 3) || (T.rows() == 3 && T.cols() == 4) && "expecting T to be 4 x 3 matrix or a 3x 4 matrix");
	if (T.rows() == 3 && T.cols() == 4)
	{
		T = T.transpose().eval();
	}
	int num_b = P0.rows() / 4;
	for (int i = 0; i < num_b; i++)
	{
		P0.block(4 * i, 0, 4, 3) = T;
	}
}

/*
Sets rig parameters  corresponding to bone BI to the affine transform T
*/
void set_rig_parameters(MatrixXd& P0, MatrixXd& T, int BI)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	assert((T.rows() == 4 && T.cols() == 3) || (T.rows() == 3 && T.cols() == 4) && "expecting T to be 4 x 3 matrix or a 3x 4 matrix");
	assert(BI < (P0.rows() / 4) && "Bone index is higher than the number of bones given!");
	if (T.rows() == 3 && T.cols() == 4)
	{
		T = T.transpose().eval();
	}
	int i = BI * 4;
	P0.block(i, 0, 4, 3) = T;
}
/*
Sets flattened parameters  corresponding to bone BI to the affine transform T
*/
void set_rig_parameters(VectorXd& p0, MatrixXd& T, int BI)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P0 = Map<MatrixXd>(p0.data(), num_b * (dim + 1), dim);
	set_rig_parameters(P0, T, BI);
	p0 = Map<VectorXd>(P0.data(), (dim * (dim + 1)) * num_b);
}

/*
Sets all rig parameters to the T affine matrix on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void set_rig_parameters_identity(VectorXd& p0, MatrixXd& T)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P0 = Map<MatrixXd>(p0.data(), num_b * (dim + 1), dim);
	set_rig_parameters(P0, T);
	p0 = Map<VectorXd>(P0.data(), (dim * (dim + 1)) * num_b);
}

/*
Sets all rig parameters to the identity affine matrix
*/
void set_rig_parameters_identity(MatrixXd& P0)
{
	assert(P0.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P0.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	int num_b = P0.rows() / 4;
	for (int i = 0; i < num_b; i++)
	{
		P0.block(4 * i, 0, 3, 3) = Matrix3d::Identity();
	}
}

/*
Sets all rig parameters to the identity affine matrix on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void set_rig_parameters_identity(VectorXd& p0)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P0 = Map<MatrixXd>(p0.data(), num_b * (dim + 1), dim);
	set_rig_parameters_identity(P0);
	p0 = Map<VectorXd>(P0.data(), (dim * (dim + 1)) * num_b);
}

void get_translations_from_rig_parameters(MatrixXd& P, MatrixXd& T)
{
	assert(P.cols() == 3 && "expecting stacked (4 #bones) x 3 P0");
	assert((P.rows() % 4) == 0 && "expecting stacked (4 #bones) x 3 P0");
	int num_b = P.rows() / 4;
	T.resize(num_b, 3);
	for (int i = 0; i < num_b; i++)
	{
		T.row(i) = P.row(4 * i + 3);
	}
}

void get_translations_from_rig_parameters(VectorXd& p, MatrixXd& T)
{
	assert((p.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	int num_b = p.rows() / 12;
	int dim = 3;
	MatrixXd P = Map<MatrixXd>(p.data(), num_b * (dim + 1), dim);
	get_translations_from_rig_parameters(P, T);
}


void world_to_rel_rig_parameters(const VectorXd& p_w, const VectorXd& p0, VectorXd& p_rel)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	assert(p0.rows() == p_w.rows());
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P_w = Map<const MatrixXd>(p_w.data(), num_b * 4, 3);
	MatrixXd P0 = Map<const MatrixXd>(p0.data(), num_b * 4, 3);
	MatrixXd P_rel = MatrixXd::Zero(P0.rows(), P0.cols());

	Matrix4d A_w = Matrix4d::Identity();
	Matrix4d A0 = Matrix4d::Identity();
	MatrixXd A_rel;
	for (int i = 0; i < num_b; i++)
	{
		A_w.block(0, 0, 3, 4) = P_w.block(i * 4, 0, 4, 3).transpose();
		A0.block(0, 0, 3, 4) = P0.block(i * 4, 0, 4, 3).transpose();
		Matrix4d Ai = A0.inverse();

		A_rel = ((A_w * Ai).block(0, 0, 3, 4)).transpose();
		P_rel.block(i * 4, 0, 4, 3) = A_rel;
	}
	p_rel = Map<VectorXd>(P_rel.data(), p0.rows());
}

void rel_to_world_rig_parameters(const VectorXd& p_rel, const VectorXd& p0, VectorXd& p_w)
{
	assert((p0.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	assert(p0.rows() == p_rel.rows());
	int num_b = p0.rows() / 12;
	int dim = 3;
	MatrixXd P_rel = Map<const MatrixXd>(p_rel.data(), num_b * 4, 3);
	MatrixXd P0 = Map<const MatrixXd>(p0.data(), num_b * 4, 3);
	MatrixXd P_w = MatrixXd::Zero(P0.rows(), P0.cols());

	Matrix4d A_rel = Matrix4d::Identity();
	Matrix4d A0 = Matrix4d::Identity();
	MatrixXd A_w;
	for (int i = 0; i < num_b; i++)
	{
		A_rel.block(0, 0, 3, 4) = P_rel.block(i * 4, 0, 4, 3).transpose();
		A0.block(0, 0, 3, 4) = P0.block(i * 4, 0, 4, 3).transpose();
	

		A_w = (A_rel * A0).block(0, 0, 3, 4).transpose();//  ((A_w).block(0, 0, 3, 4)).transpose();
		P_w.block(i * 4, 0, 4, 3) = A_w;
	}
	p_w = Map<VectorXd>(P_w.data(), p0.rows());
}

/*
Returns 3x4 affine bone matrix for bone i from flattened rig parameters p
*/
void get_bone_transform(const VectorXd& p, int i, MatrixXd& A)
{
	assert((p.rows() % (12)) == 0 && "flattened rig parameters in 3D mush have p0.rows() be a multiple of 12");
	A.resize(3, 4);
	A.setZero();
	MatrixXd P = Map<const MatrixXd>(p.data(), p.rows() / 3, 3);
	A = P.block(4 * i, 0, 4, 3).transpose();
}
MatrixXd get_bone_transform(const VectorXd& p, int i)
{
	MatrixXd A;
	get_bone_transform(p, i, A);
	return A;
}

void world_to_rel_rig_anim(MatrixXd& P_w, VectorXd& p0, MatrixXd& P_r)
{
	assert(P_w.rows() == p0.rows() && "Rig anim has different #parameters than p0");
	assert(P_w.cols() > 0 && "need at least one frame of animation to covnert");
	int num_b = P_w.rows() / 12;
	int dim = 3;

	int nf = P_w.cols();
	VectorXd p_w, p_rel;
	P_r.resize(P_w.rows(), P_w.cols());
	P_r.setZero();
	for (int i = 0; i < nf; i++)
	{
		p_w = P_w.col(i);
		world_to_rel_rig_parameters(p_w, p0, p_rel);
		P_r.col(i) = p_rel;
	}

}

