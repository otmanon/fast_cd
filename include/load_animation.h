#pragma once
#include <Eigen/Core>
/*
Loads an animation file found in the anim_filepath.json file. 
Flattens all rig parameters (row-wise), and coincatenates them column wise into anim_P
*/
void load_animation(std::string anim_filepath, Eigen::MatrixXd& anim_P, bool is_global);


/*
Same as above, but also attempts to fit the global transformations, scale, rotations to a rest p0, with parent pI and bone length
*/
void load_animation_and_fit_to_rest_skeleton_json(std::string anim_filepath, Eigen::VectorXd& p0, Eigen::VectorXi& pI, Eigen::MatrixXd& anim_P, bool& is_global);


void load_animation_json(std::string anim_filepath, Eigen::MatrixXd& anim_P, bool& is_global);