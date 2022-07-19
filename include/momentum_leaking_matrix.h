#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace fast_cd
{

    enum MOMENTUM_LEAK_MATRIX
    {
        MOMENTUM_LEAK_DIFFUSION= 0,
        MOMENTUM_LEAK_DIFFUSION_FLIP = 1,
    };
}

using namespace Eigen;
/*
Construct the momentum leaking matrix D for complementary dynamics. For normal CD, momentum wouldn't leak from the rig to the mesh,
so we introduce D that softens our complementarity enforcement. D is obtained through a simple diffusion of the shape, such that
vertices closer to the surface leak more momentum, or less momentum, than those in the interior.
*/
void momentum_leaking_matrix(const MatrixXd& V, const MatrixXi& F, fast_cd::MOMENTUM_LEAK_MATRIX method, SparseMatrix<double>& D, double dt = 1e-6);