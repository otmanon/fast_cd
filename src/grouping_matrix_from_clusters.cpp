#include "grouping_matrix_from_clusters.h"
#include "interweaving_matrix.h"
#include <igl/colon.h>
#include <igl/repdiag.h>
void grouping_matrix_from_clusters(Eigen::VectorXi& I, Eigen::SparseMatrix<double>& G)
{
    Eigen::VectorXi J;
    igl::colon(0, I.rows()-1, J);

    typedef Eigen::Triplet<double> T;

    std::vector<T> tripletList(I.rows());
    for (int i = 0; i < I.rows(); i++)
    {
        tripletList.emplace_back(T(I(i), J(i), 1));
    }
    G.resize(I.maxCoeff() + 1.0, J.rows());
    G.setFromTriplets(tripletList.begin(), tripletList.end());

}