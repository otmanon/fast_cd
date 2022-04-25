#include "kmeans_jacobi.h"
#include <igl/AABB.h>
#include <igl/LinSpaced.h>

// Jacobi-style online updates. These are worst-case costly because a single
// iteration through all points could be O(n k log k).
// 
// Inputs:
//   P  #P by dim list of input points
//   W  #P list of positive weights
//   k  number of clusters
//   C  k by dim list of cluster centroids initial guesses
//   I  #P list of indices into rows of C
//   sqrD  #P list of squared distances to corresponding cluster center
//   Wcount  k list of weighted counts
//   obj  current scalar objective function 
// Outputs:
//   C,I,sqr,Wcount,obj  updated to converged values
void igl::kmeans_jacobi(
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & W,
  const int k,
  Eigen::MatrixXd & C,
  Eigen::VectorXi & I,
  Eigen::VectorXd & sqrD,
  Eigen::VectorXd & Wcount,
  double & obj)
{
  const int n = P.rows();
  assert(C.rows() == k && "C should contain initial guesses");
  assert(I.size() == P.rows() && "I should contain initial guesses"); 
  assert(sqrD.size() == P.rows() && "I should contain initial guesses"); 
  // obj, sqrD, and I have just been updated.
  // sigh. AABB is templated on dim directly, so we have to duplicate code here
  // (or propagate dim template to this function). This will likely fail to
  // compile if ColsAtCompileTime for P is ≠ Dynamic,2,3.
  igl::AABB<Eigen::MatrixXd,3> tree3;
  igl::AABB<Eigen::MatrixXd,2> tree2;
  // for AABB tree (if used)
  Eigen::VectorXi Ele = igl::LinSpaced<Eigen::VectorXi>(k,0,k-1);
  switch(P.cols())
  {
    default: break;
    case 3:
      tree3.init(C,Ele);
      break;
    case 2:
      tree2.init(C,Ele);
      break;
  }

  // Q: If the batch iterations "converged" then are these necessary? (cost of a
  // single iteration to confirm is not high in any case)
  const int max_individual_iter = 100;
  int individual_iter = 0;
  while(true)
  {
    //printf("O: obj(%d): %g\n",individual_iter,obj);
    const double prev_individual_iter_obj = obj;
    bool any_changed = false;
    // consider each point O(n k log k) because of possible tree rebuild...
    for(int j = 0;j<n;j++)
    {
      // find closest cluster center to point j
      const int prev_Ij = I(j);
      const double prev_sqrDj = sqrD(j);
      // O(k) or O(log k)
      switch(P.cols())
      {
        default: 
        {
          for(int i = 0;i<k;i++)
          {
            const double sqrdji = (C.row(i)-P.row(j)).squaredNorm();
            if(sqrdji < sqrD(j))
            {
              sqrD(j) = sqrdji;
              I(j) = i;
            }
          }
            break;
        }
        case 2:
        {
          const Eigen::RowVector2d Pj = P.row(j);
          Eigen::RowVector2d _;
          sqrD(j) = tree2.squared_distance(C,Ele,Pj,I(j),_);
          break;
        }
        case 3:
        {
          const Eigen::RowVector3d Pj = P.row(j);
          Eigen::RowVector3d _;
          sqrD(j) = tree3.squared_distance(C,Ele,Pj,I(j),_);
          break;
        }
      }
      // should we switch point j?
      if(prev_Ij != I(j))
      {
        const auto WPj = (W(j)*P.row(j)).eval();
        // downdate
        C.row(prev_Ij) = (C.row(prev_Ij)*Wcount(prev_Ij))-WPj;
        Wcount(prev_Ij) -= W(j);
        C.row(prev_Ij) /= Wcount(prev_Ij);
        // udpate
        C.row(I(j)) = (C.row(I(j))*Wcount(I(j)))+WPj;
        Wcount(I(j)) += W(j);
        C.row(I(j)) /= Wcount(I(j));
        obj = obj + W(j)*sqrD(j) - W(j)*prev_sqrDj;
        any_changed = true;
        // tree is now out of date... rebuild. O(k log k)
        switch(P.cols())
        {
          default: break;
          case 3: tree3.init(C,Ele); break;
          case 2: tree2.init(C,Ele); break;
        }
      }
    }
    if(prev_individual_iter_obj==obj)
    {
      //printf("converged @ %d : %g\n",individual_iter,obj);
      break;
    }else if(prev_individual_iter_obj>obj)
    {
      // this should not happen.
      assert(false);
    }else if(individual_iter+1 == max_individual_iter)
    {
      //printf("failed to converge @ %d\n", individual_iter);
      break;
    }
    individual_iter++;
  }
  printf("[kmeans] individual_iter: %d\n",individual_iter);
}
