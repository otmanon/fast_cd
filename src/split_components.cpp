#include "split_components.h"
#include <queue>

void igl::split_components(
  const Eigen::Array<bool,Eigen::Dynamic,1> & I,
  const Eigen::SparseMatrix<bool> & A,
  std::vector<std::vector<int> > & S)
{
  typedef Eigen::Index Index;
  const auto m = A.rows();
  assert(A.cols() == A.rows() && "A should be square");
  std::vector<bool> seen(m);
  // O(n)
  for(int i = 0;i<I.size();i++)
  {
    bool new_list = true;
    std::queue<int> Q;
    Q.push(i);
    while(!Q.empty())
    {
      const int g = Q.front();
      Q.pop();
      // already seen
      if(seen[g]) continue;
      // dead
      if(!I(g)) continue;
      // see it
      seen[g] = true;
      if(new_list)
      {
        S.push_back({});
        new_list = false;
      }
      S.back().push_back(g);
      for(Eigen::SparseMatrix<bool>::InnerIterator it (A,g); it; ++it)
      {
        const Index n = it.index();
        // already seen
        if(seen[n]) continue;
        // dead
        if(!I(n)) continue;
        Q.push(n);
      }
    }
  }
}

void igl::split_components(
  const Eigen::VectorXi & I,
  const Eigen::SparseMatrix<bool> & A,
  Eigen::VectorXi & J)
{
  int k = 0;
  const int maxI = I.maxCoeff();
  J.setZero(I.size());
  for(int i =0;i<=maxI;i++)
  {
    std::vector<std::vector<int>> S;
    split_components( (I.array()==i).eval(), A, S);
    for(auto & s : S)
    {
      for(auto & j : s)
      {
        J(j) = k;
      }
      k++;
    }
  }
}

#include <igl/adjacency_matrix.h>

void igl::split_components(
  const Eigen::VectorXi & I,
  const Eigen::MatrixXi & F,
  Eigen::VectorXi & J)
{
  Eigen::SparseMatrix<bool> A;
  igl::adjacency_matrix(F,A);
  return split_components(I,A,J);
}


void igl::split_components(
  const Eigen::MatrixXd & B,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & C)
{
  Eigen::SparseMatrix<bool> A;
  igl::adjacency_matrix(F,A);
  int k = 0;
  C.resize(B.rows(),B.cols());
  const int n = B.cols();
  for(int i =0;i<n;i++)
  {
    std::vector<std::vector<int>> S;
    split_components( (B.col(i).array()>0).eval(), A, S);
    // make room
    C.conservativeResize(C.rows(),C.cols()*2+S.size());
    for(auto & s : S)
    {
      C.col(k).setZero();
      for(auto & j : s)
      {
        C(j,k) = B(j,i);
      }
      k++;
    }
  }
  C.conservativeResize(C.rows(),k);
}
