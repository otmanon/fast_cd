#pragma once
#include <Eigen/Core>
#include <filesystem>
#include <igl/writeDMAT.h>

struct sim_metrics
{
    Eigen::VectorXd inertias, volumes, bendings;
    Eigen::VectorXd sim_timings, proj_timings;
    void init(int n)
    {
        inertias.resize(n);
        volumes.resize(n);
        bendings.resize(n);
        sim_timings.resize(n);
        proj_timings.resize(n);
    }
    void update(int i, double inertia, double volume, double bending, double sim_timing, double proj_timing)
    {
        inertias(i) = inertia;
        volumes(i) = volume;
        bendings(i) = bending;
        sim_timings(i) = sim_timing;
        proj_timings(i) = proj_timing;
    }

    void write(std::string results_path)
    {
        if (!std::filesystem::exists(std::filesystem::path(results_path).parent_path()))
        {
            std::filesystem::create_directories(std::filesystem::path(results_path).parent_path());
        }
        Eigen::MatrixXd results(inertias.rows(), 5);
        results.col(0) = inertias; results.col(1) = volumes; results.col(2) = bendings; results.col(3) = sim_timings; results.col(4) = proj_timings;
        igl::writeDMAT(results_path, results);
    }

} m_reduced, m_full;
