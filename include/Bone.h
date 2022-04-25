#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <igl/PI.h>
class Bone
{
    const double inf = std::numeric_limits<double>::infinity();

public:
    //Constructor
    Bone(
        int _parent_index,
        int _weight_index,
        Eigen::Affine3d _rest_T,
        double _length, Eigen::VectorXi _wI, Eigen::VectorXd _wV);
public:
    // Index into skeleton bone-list  of parent (-1 indicates root)
    int parent_index;
    // Index into columns of weights matrix for corresponding linear blend
    // skinning weights (-1 indicates no associated weights)
    int weight_index;
    // Transformation mapping from "canonical bone" to "rest bone"
    Eigen::Affine3d rest_T;
    // Length of bone.
    double length;
    // Euler Angles of current pose  // this is currently not used
    Eigen::Vector3d xzx;
    // max, min Euler Angles: joint limits
    Eigen::Vector3d xzx_max, xzx_min;


    //why not just remember this here!
    //weight Indeces, the indeces into V that this bone affects.
    Eigen::VectorXi wI;
    //weight values. The amount this bone affects the corresponding vertex wI
    Eigen::VectorXd wV;

};
