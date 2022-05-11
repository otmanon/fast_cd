#pragma once
#include <Eigen/Core>

struct sim_params {
    Eigen::VectorXd z_curr, z_prev;
    Eigen::VectorXd p_curr, p_prev;

    void init(Eigen::VectorXd& z, Eigen::VectorXd& p)
    {
        this->z_curr = z;
        this->z_prev = z;

        this->p_curr = p;
        this->p_prev = p;
    }

    void update(Eigen::VectorXd& z, Eigen::VectorXd& p)
    {
        this->p_prev = p_curr;
        this->p_curr = p;

        this->z_prev = z_curr;
        this->z_curr = z;
    }
} s;