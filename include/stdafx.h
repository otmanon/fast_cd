#pragma once
#include <thread>

//common igl functions
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/project.h>
#include "igl/repdiag.h"
#include "igl/cat.h"
#include "igl/polar_svd3x3.h"
#include "igl/cotmatrix.h"
#include "igl/arap_linear_block.h"
#include "igl/covariance_scatter_matrix.h"
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice_into.h>
#include "igl/repmat.h"
#include "igl/colon.h"


#include <stdio.h>
#include <assert.h>
#include <iostream>

#include "Spectra/MatOp/SymShiftInvert.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Spectra/SymGEigsShiftSolver.h"
#include "Spectra/SymGEigsSolver.h"
#include "Spectra/MatOp/SparseRegularInverse.h"



#include <Eigen/sparse>
#include <Eigen/core>