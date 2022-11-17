#pragma once
#include <Eigen/Core>
using namespace Eigen;
/*
Translates _stacked_ affine rig parameters in P0 all by t
*/
void translate_rig_parameters(MatrixXd& P0, Vector3d& t);


/*
Transforms _stacked_ affine rig parameters in P0 all by T
*/
void transform_rig_parameters(MatrixXd& P0, MatrixXd& T);

void transform_rig_parameters(VectorXd& p, MatrixXd& A);

/*
Transforms a rig animation. Given a stacked list of flattened rig parameters, apply the transform A to them
*/
void transform_rig_parameters_anim(MatrixXd& anim_P, MatrixXd& A);

/*
Translates _stacked_ affine rig parameters in P0 all by t on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void rotate_rig_parameters(VectorXd& p0, Vector3d& t);

/*
Rotates _stacked_ affine rig parameters in P0 all by R
*/
void rotate_rig_parameters(MatrixXd& P0, Matrix3d& R);

/*
Rotates _stacked_ affine rig parameters in P0 all by R on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void rotate_rig_parameters(VectorXd& p0, Matrix3d& R);

/*
Sets all rig parameters to the same affine transform T
*/
void set_rig_parameters(MatrixXd& P0, MatrixXd& T);

/*
Sets rig parameters  corresponding to bone BI to the affine transform T
*/
void set_rig_parameters(MatrixXd& P0, MatrixXd& T, int BI);
/*
Sets flattened parameters  corresponding to bone BI to the affine transform T
*/
void set_rig_parameters(VectorXd& p0, MatrixXd& T, int BI);
/*
Sets all rig parameters to the T affine matrix on flattened input/outpout.
Assumes columnb order flattening convenction
*/
void set_rig_parameters_identity(VectorXd& p0, MatrixXd& T);

/*
Sets all rig parameters to the identity affine matrix
*/
void set_rig_parameters_identity(MatrixXd& P0);

/*
Sets all rig parameters to the identity affine matrix on flattened input/outpout. 
Assumes columnb order flattening convenction 
*/
void set_rig_parameters_identity(VectorXd& p0);

void get_translations_from_rig_parameters(MatrixXd& P, MatrixXd& T);

void get_translations_from_rig_parameters(VectorXd& p, MatrixXd& T);


void world_to_rel_rig_parameters(const VectorXd& p_w, const VectorXd& p0, VectorXd& p_rel);

void rel_to_world_rig_parameters(const VectorXd& p_rel, const VectorXd& p0, VectorXd& p_w);
/*
Returns 3x4 affine bone matrix for bone i from flattened rig parameters p
*/
void get_bone_transform(const VectorXd& p, int i, MatrixXd& A);
MatrixXd get_bone_transform(const VectorXd& p, int i);

void world_to_rel_rig_anim(MatrixXd& P_w, VectorXd& p0, MatrixXd& P_r);

