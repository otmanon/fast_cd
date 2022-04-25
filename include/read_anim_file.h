#pragma once
#include <Skeleton.h>
#include <filesystem>
#include <string>

std::vector<std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> > read_anim_file(const std::string file, const Skeleton& skeleton , Eigen::Matrix4d& S, Eigen::Vector3d& t);