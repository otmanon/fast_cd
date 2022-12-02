#pragma once

#include "get_modes.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <iostream>
#include <filesystem>
using namespace Eigen;

/*
Reads/Writes cached modes by looking in modes_dir. If B.DMAT or L.DMAT can't be found, both are recomputed
Only looks in cache if cache flag is set to true. 

Input:
V -> n x 3 geometry
T -> f x 4 tets
W -> n x b weights
J -> 3n x 12b lbs jacobian, may act as constraint on modes
mode_type -> string specifiying mode type ("displacement" or "skinning")
modes_dir -> string, directory where we will search for modes and write to. If it doesn't exist, will be created
cache -> boolean, set to true if you want to load/write modes. otherwise will compute from scratch no matter what
num_modes -> int, number of modes to search for/compute. If cache have less modes than this, will recompute

Output
B -> either 3n x num_modes for displacement modes, or n x num_modes/12 for skinning modes
L -> either num_modes vector or num_modes/12 vector of eigenvalues
Ws -> if we picked skinning subspace, returns Ws
*/
void get_cached_modes(MatrixXd& V, MatrixXi& T, MatrixXd& W, SparseMatrix<double>& J, std::string mode_type, 
    std::string modes_dir, bool cache, int num_modes, MatrixXd& B, VectorXd& L, MatrixXd& Ws)
{
    namespace fs = std::filesystem;
    B.resize(0, 0);
    L.resize(0);
    std::string dir = modes_dir + "/" + mode_type + "/";
    bool found = false;
    if (cache)
    {
        std::cout << "looking for cached modes in " << dir << std::endl;

        found = igl::readDMAT(dir + "B.DMAT", B) && igl::readDMAT(dir + "L.DMAT", L);
        if (mode_type == "skinning")
            found = found && igl::readDMAT(dir + "W.DMAT", Ws);
    
        if (B.cols() < num_modes || L.rows() < num_modes)
        {
            found = false;
            std::cout << "could not find adequate modes in " << dir << std::endl;
        }
        if (found)
            std::cout << "found them !" << std::endl;
    }
    if (!found)
    {
        std::cout << "recomputing " << num_modes << " modes from scratch ..." << std::endl;
        get_modes(V, T, J, mode_type, num_modes, B, L, Ws);
      
        std::cout << "saving modes in " << dir <<  "..." << std::endl;

            if (!fs::exists(fs::path(dir)))
                fs::create_directories(fs::path(dir));
            igl::writeDMAT(dir + "B.DMAT", B); igl::writeDMAT(dir + "L.DMAT", L);
            igl::writeDMAT(dir + "W.DMAT", Ws);
      
    }
}