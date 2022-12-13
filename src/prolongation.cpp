#include "prolongation.h"

#include <igl/point_mesh_squared_distance.h>
#include <igl/barycentric_coordinates.h>
#include <igl/slice.h>
#include <igl/boundary_facets.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/in_element.h>
#include <igl/AABB.h>

void prolongation(const Eigen::MatrixXd& X_fine, const Eigen::MatrixXd& X_coarse, const Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W)
{
	Eigen::VectorXi I;

	igl::AABB<Eigen::MatrixXd, 3> tree;
	tree.init(X_coarse, T_coarse);
	//igl::point_mesh_squared_distance(X_fine, X_coarse, T_coarse, sqrD, I, CP);

	igl::in_element(X_coarse, T_coarse, X_fine, tree, I);
	Eigen::MatrixXi T_fine;
	
	Eigen::MatrixXi F;
	Eigen::VectorXi n, FiT;
	igl::boundary_facets(T_coarse, F,FiT, n);

	Eigen::MatrixXd CP;
	CP = X_fine;			//will only be not true for fine elements poking out of coarse mesh
	int ind, ind_face; Eigen::RowVector3d cp;
	tree.init(X_coarse, F);				//need to reinitialize this 
	int ci = 0;
	for (int i = 0; i < I.rows(); i++)
	{
		
		if (I(i) < 0)

		{ //need to project to closest point on mesh surface
			tree.squared_distance(X_coarse, F, X_fine.row(i), ind_face, cp); //convert this to surface mesh
			ind = FiT(ind_face);
			CP.row(i) = cp;
			I(i) = ind;
			ci++;
		}
	}
	
	
	
	igl::slice(T_coarse, I, 1, T_fine);  //slice out coarse tets for each fine vertex

	Eigen::VectorXi A, B, C, D;			//indices of first, second, third and fourth vertex in each tet.
	Eigen::MatrixXd XA, XB, XC, XD;    // vertex positions for each of the above vertices
	A = T_fine.col(0); igl::slice(X_coarse, A, 1, XA);
	B = T_fine.col(1); igl::slice(X_coarse, B, 1, XB);
	C = T_fine.col(2); igl::slice(X_coarse, C, 1, XC);
	D = T_fine.col(3); igl::slice(X_coarse, D, 1, XD);

	Eigen::MatrixXd BC; //barycentric coordinates... can use these barycentric coordinates to make an embedding matrix W
	igl::barycentric_coordinates(CP, XA, XB, XC, XD, BC);

	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(X_fine.rows() * 4);
	for (int i = 0; i < X_fine.rows(); i++)
	{
		tripletList.emplace_back(i, A(i), BC(i, 0));
		tripletList.emplace_back(i, B(i), BC(i, 1));
		tripletList.emplace_back(i, C(i), BC(i, 2));
		tripletList.emplace_back(i, D(i), BC(i, 3));
	}

	//std::cout << BC << std::endl;
	W.resize(X_fine.rows(), X_coarse.rows());
	W.setFromTriplets(tripletList.begin(), tripletList.end());
}



void prolongation(const Eigen::MatrixXd& X_fine,const  Eigen::MatrixXd& X_coarse,const  Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W, Eigen::SparseMatrix<double>& WI)
{
	Eigen::VectorXi I;

	igl::AABB<Eigen::MatrixXd, 3> tree;
	tree.init(X_coarse, T_coarse);
	//igl::point_mesh_squared_distance(X_fine, X_coarse, T_coarse, sqrD, I, CP);

	igl::in_element(X_coarse, T_coarse, X_fine, tree, I);
	Eigen::MatrixXi T_fine;

	Eigen::MatrixXi F;
	Eigen::VectorXi n, FiT;
	igl::boundary_facets(T_coarse, F, FiT, n);

	Eigen::MatrixXd CP;
	CP = X_fine;			//will only be not true for fine elements poking out of coarse mesh
	int ind, ind_face; Eigen::RowVector3d cp;
	tree.init(X_coarse, F);				//need to reinitialize this 
	int ci = 0;
	for (int i = 0; i < I.rows(); i++)
	{

		if (I(i) < 0)

		{ //need to project to closest point on mesh surface
			tree.squared_distance(X_coarse, F, X_fine.row(i), ind_face, cp); //convert this to surface mesh
			ind = FiT(ind_face);
			CP.row(i) = cp;
			I(i) = ind;
			ci++;
		}
	}
	igl::slice(T_coarse, I, 1, T_fine);  //slice out coarse tets for each fine vertex

	Eigen::VectorXi A, B, C, D;			//indices of first, second, third and fourth vertex in each tet.
	Eigen::MatrixXd XA, XB, XC, XD;    // vertex positions for each of the above vertices
	A = T_fine.col(0); igl::slice(X_coarse, A, 1, XA);
	B = T_fine.col(1); igl::slice(X_coarse, B, 1, XB);
	C = T_fine.col(2); igl::slice(X_coarse, C, 1, XC);
	D = T_fine.col(3); igl::slice(X_coarse, D, 1, XD);

	Eigen::MatrixXd BC; //barycentric coordinates... can use these barycentric coordinates to make an embedding matrix W
	igl::barycentric_coordinates(CP, XA, XB, XC, XD, BC);

	std::vector<Eigen::Triplet<double>> tripletList, tripletListI;
	tripletList.reserve(X_fine.rows() * 4);
	for (int i = 0; i < X_fine.rows(); i++)
	{
		tripletList.emplace_back(i, A(i), BC(i, 0));
		tripletList.emplace_back(i, B(i), BC(i, 1));
		tripletList.emplace_back(i, C(i), BC(i, 2));
		tripletList.emplace_back(i, D(i), BC(i, 3));

		tripletListI.emplace_back(i, A(i), 1.0);
		tripletListI.emplace_back(i, B(i), 1.0);
		tripletListI.emplace_back(i, C(i), 1.0);
		tripletListI.emplace_back(i, D(i), 1.0);
	}

	//std::cout << BC << std::endl;
	W.resize(X_fine.rows(), X_coarse.rows());
	W.setFromTriplets(tripletList.begin(), tripletList.end());

	WI.resize(W.rows(), W.cols());//this matrix contains the same stuff as W, but only has 1's on its nonzero terms
	WI.setFromTriplets(tripletListI.begin(), tripletListI.end());

	//Eigen::MatrixXd diff = X_fine - W * X_coarse;
}