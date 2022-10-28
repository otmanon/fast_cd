#pragma once
#include <Eigen/Core>

using namespace Eigen;
/*
Taking the LBS product C'B(W), where B is the lbs jacobian, we can rewrite it a reshaped product directly on W. When C itself 
has the structure of an lbs matrix, a lot of redundant rows appear. Luckily they appear in a very predictable way. 

This returns the non-redundant rows of the matrix A, in AW, 
Input
num_b -> number of bones giving rise to the lbs matrix C.
dim -> dimensionality of the problem (2 or 3)

Output
VectorXi -> indices into A of non-redundant rows

*/
void redundant_indices_JV(int num_b, int dim, VectorXi& I)
{
	I.resize(0);
	for (int oi = 1; oi <= dim + 1; oi++)
	{
		int offset = num_b * dim * (dim + 1) * (oi - 1);
			for (int b = 1; b <= num_b; b++)
			{
				offset += (dim + 1) * (b - 1);
				for (int i = oi; i <= dim + 1; i++)
				{
					I.conservativeResize(I.rows() + 1);
					I(I.rows() - 1) = i + offset;
				}
			}
	}

	I = I.array() - 1;
	// this follows the matlab syntax, so just remove 1 from index at the end
	
}
