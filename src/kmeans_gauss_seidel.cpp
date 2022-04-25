#include "kmeans_gauss_seidel.h"
#include "igl/unique.h"
bool igl::kmeans_gauss_seidel(
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & W,
  const int k,
  const std::function<void(const Eigen::MatrixXd&,const Eigen::MatrixXd&,Eigen::VectorXi&,Eigen::VectorXd&)> & nearest_neighbor,
  Eigen::MatrixXd & C,
  Eigen::VectorXi & I,
  Eigen::VectorXd & sqrD,
  Eigen::VectorXd & Wcount,
  double & obj)
{
  assert(P.rows() == W.rows());
  assert(P.rows() >= k);
  // number of data points
  const int n = P.rows();
  // batch updates (it is known that this will not necessary converge to a local
  // minimum, but are these guaranteed at least to converge?) If yes, then could
  // replace with while loop.
  //
  // I no longer understand the statement/question above. If it converges in the
  // sense of `I` no longer changing, then this must be a local minimum, no?
  // Therefore, if "it is known that this will not necessarily converge", then
  // this must imply that it doesn't converge to a local minimum. Would be great
  // to have an example where this fails and Jacobi-polish is necessary.
  //
  // more batch iterations tends to eventually find a better local minimum
  // value: in contrast to individual iterations (phase after batch iterations)
  // converges very quickly to a (possibly poor) local minimum.
  const int max_batch_iter = 1000;
  Eigen::VectorXi count = decltype(count)::Zero(k,1);
  // initialize to bogus values
  I.setConstant(n,1,-1);
  obj = std::numeric_limits<double>::infinity();
  int batch_iter = 0;
  bool converged = false;
  // Gauss-Seidel style batch updates
  while(true)
  {
    Eigen::VectorXi prev_I = I;
    // given new C, update Ican and sqrD
    Eigen::VectorXi Ican;
    nearest_neighbor(C,P,Ican,sqrD);
    // Given Ican update I and count. O(n)
    // Can't trivially apply parallel_for here because of count. Alternative
    // would be to keep a representative for each cluster and "undo" changing
    // representative if a cluster ever dips to zero.
    for(int j = 0;j<n;j++)
    {
      if((I(j)==-1) || count(I(j)>1))
      {
        if(I(j)>-1)
        {
          count(I(j))--;
        }
        I(j) = Ican(j);
        count(I(j))++;
      }
    }
    if(I==prev_I){ converged = true;break; }

#ifndef NDEBUG
    {
      Eigen::VectorXi U;
      igl::unique(I,U);
      assert(k == U.rows());
    }
#endif

    // Given I, update C
    C.setZero(k,C.cols());
    Wcount.setZero(k);
    for(int j = 0;j<n;j++)
    {
      C.row(I(j)) += W(j)*P.row(j);
      Wcount(I(j)) += W(j);
    }
    for(int i = 0;i<k;i++)
    {
      C.row(i) /= Wcount(i);
    }

    // These stopping criteria are rarely (some never) hit.
    const double prev_obj = obj;
    obj = W.transpose() * sqrD;
    if(obj==prev_obj)
    {
      break;
    }else if(obj>=prev_obj)
    {
      break;
    }
    else if(batch_iter+1 == max_batch_iter)
    {
      break;
    }
    batch_iter++;
  }
  //printf("[kmeans] batch_iter: %d\n",batch_iter);
  return converged;
}
