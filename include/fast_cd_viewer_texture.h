#pragma once
#include "fast_cd_viewer.h"
#include "prolongation.h"

#include <igl/boundary_facets.h>
#include <igl/readOBJ.h>
#include <igl/repdiag.h>
using namespace std;
using namespace igl;
struct fast_cd_viewer_texture : public fast_cd_viewer
{

	bool texture_on;

	SparseMatrix<double> P;

	MatrixXd Xc; MatrixXi Fc; //coarse mesh
	MatrixXd Xf; MatrixXi Ff; //fine mesh

	fast_cd_viewer_texture() : fast_cd_viewer()
	{
		texture_on = false;
	};

	void set_mesh(MatrixXd& V, MatrixXi& F, int id)
	{
		Xc = F; Fc = F;
		fast_cd_viewer::add_mesh(V, F, id);
	}
	/*
	Sets the viewer to display mesh V with connectivity F, but with texture info specified
	in texture_obj and texture_img. 
	*/
	void set_mesh(MatrixXd& V, MatrixXi& T, string texture_obj, string texture_png,
		 int id, double so=1.0, RowVector3d to=RowVector3d(0, 0, 0))
	{
		Xc = V; 
		igl::boundary_facets(T, Fc);

		MatrixXd  N; MatrixXi FN;
		MatrixXd TC; MatrixXi FTC;
		readOBJ(texture_obj, Xf, TC, N, Ff, FTC, FN);
		
		Xf = (Xf * so).rowwise() - to;
		fast_cd_viewer::set_mesh(Xf, Ff, id);
		fast_cd_viewer::set_texture(texture_png, TC, FTC, id);
		prolongation(Xf, V, T, P);
		P = igl::repdiag(P, 3);
		
		texture_on = true;
		//fast_cd_viewer::invert_normals(, 0);
	}

	/*
	Displays fine mesh with coarse displacement field U.
	If viewer is not configured with a texture, sets it to coarse mesh.
	*/
	void set_coarse_vertices(MatrixXd& U, bool id=0)
	{
		MatrixXd V;
		if (texture_on)
		{
			VectorXd u = Map<VectorXd>(U.data(),
				U.rows()*U.cols(), 1);

			VectorXd v = P * u;
			V = Map<MatrixXd>(v.data(), v.rows() / 3, 3) + Xf;	
		}
		else
		{
			V = U + Xc;
		}
		fast_cd_viewer::set_vertices(V, id);
	}

	/*
Displays fine mesh with coarse displacement field U.
If viewer is not configured with a texture, sets it to coarse mesh.
*/
	void set_coarse_vertices(VectorXd& u, bool id = 0)
	{
		MatrixXd V;
		if (texture_on)
		{
			VectorXd v = P * u;
			V = Map<MatrixXd>(v.data(), int(v.rows() / 3), 3) + Xf;
		}
		else
		{
			MatrixXd U = Map<MatrixXd>(u.data(), int(u.rows()/3),3);;
			V = U + Xc;
		}
		fast_cd_viewer::set_vertices(V, id);
	}

	/*
	Sets the viewer to display mesh with connectivity F,
	but with texture info specified in texture_obj and texture_png
	*/
	//void set_mesh(MatrixXd& V, MatrixXi&  F, MatrixXd& Vf, MatrixXi& Ff, MatrixXd& TC, MatrixXi& FTC)
	//{

	//}

};