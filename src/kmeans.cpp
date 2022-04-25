#include "kmeans.h"
#include "kmeans_gauss_seidel.h"
#include "kmeans_jacobi.h"
#include <igl/unique_rows.h>
#include <igl/parallel_for.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/LinSpaced.h>
#include <igl/randperm.h>
#include <igl/get_seconds.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <igl/cumsum.h>
#include <limits>
#include <cstdio>
#include <cassert>
#include <functional>
#include <iostream>
void igl::kmeans(
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & W,
  const int k,
  const std::function<void(const Eigen::MatrixXd&,const Eigen::MatrixXd&,Eigen::VectorXi&,Eigen::VectorXd&)> & nearest_neighbor,
  Eigen::MatrixXd & C,
  Eigen::VectorXi & I)
{
  // Reduce to unique case
  Eigen::MatrixXd U;
  Eigen::VectorXi _,J;
  igl::unique_rows(P,U,_,J);
  Eigen::VectorXi IU;
  Eigen::VectorXd WU = Eigen::VectorXd::Zero(U.rows(),1);
  for(int i = 0;i<P.rows();i++)
  {
    WU(J(i))+=W(i);
  }

  const auto plus_plus_init = [&nearest_neighbor](
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C)
  {
    const int n = P.rows();
    // 1. Choose one center uniformly at random among the data points.
    int i = 0;
    Eigen::VectorXi J(k);
    J(i++) = rand() % n;
    Eigen::VectorXi I(n);
    Eigen::VectorXd sqrD = 
      Eigen::VectorXd::Constant(n,1,std::numeric_limits<double>::infinity());
    // The wikipedia "steps" imply a lot of recomputation of distances. Just
    // need to update minimum as we go.
    while(i<k)
    {
      // Update sqrD
      {
        Eigen::Matrix<double,1,Eigen::Dynamic> Clast = P.row(J(i-1));
        igl::parallel_for(n,[&](const int p)
        {
          const double new_sqrd = (P.row(p)-Clast).squaredNorm();
          if(new_sqrd < sqrD(p))
          {
            I(p) = i-1;
            sqrD(p) = new_sqrd;
          }
        },10000);
      }
      // 3. Choose one new data point at random as a new center, using a
      // weighted probability distribution where a point x is chosen with
      // probability proportional to D(x)².
      Eigen::VectorXd sqrD0 = 
        (Eigen::VectorXd(sqrD.size()+1)<<0, sqrD.array() * W.array() ).finished();
      Eigen::VectorXd CDF;
      igl::cumsum(sqrD0,1,CDF);
      const double threshold = (double)rand() / RAND_MAX * CDF(CDF.size()-1);
      int j = 0;
      for(;j<n;j++){ if(CDF(j+1)>threshold){ break; } }
      assert(j<n);
      J(i++) = j;
      // 4. Repeat Steps 2 and 3 until k centers have been chosen.
    }
    igl::slice(P,J,1,C);
  };

  const auto random_init = [](
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C)
  {
    // (uniform) random initialization → C
    Eigen::VectorXi J;
    igl::randperm(P.rows(),J);
    igl::slice(P,J.head(k).eval(),1,C);
  };

  const auto kmeans_unique = [&nearest_neighbor,&random_init,&plus_plus_init](
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I)
  {
    const auto & tictoc = []()
    {
      static double t_start = igl::get_seconds();
      double diff = igl::get_seconds()-t_start;
      t_start += diff;
      return diff;
    };
    {
      // plus_plus_init actually computes valid I,sqrD, could shave off one
      // iteration if kmeans_gauss_seidel assumed these were valid, too
      Eigen::VectorXd sqrD,Wcount;
      double obj;
      plus_plus_init(P,W,k,C);
      if(!kmeans_gauss_seidel(P,W,k,nearest_neighbor,C,I,sqrD,Wcount,obj))
      {
        kmeans_jacobi(P,W,k,C,I,sqrD,Wcount,obj);
      }
    }
  };
  kmeans_unique(U,WU,k,C,IU);

  igl::slice(IU,J,1,I);
}

void igl::kmeans(
  const Eigen::MatrixXd & P,
  const int k,
  Eigen::MatrixXd & C,
  Eigen::VectorXi & I)
{
  // uniform weights
  const Eigen::VectorXd W = Eigen::VectorXd::Constant(P.rows(),1,1);
  return kmeans(P,W,k,C,I);
}


void igl::kmeans(
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & W,
  const int k,
  Eigen::MatrixXd & C,
  Eigen::VectorXi & I)
{
  std::function<void(const Eigen::MatrixXd&,const Eigen::MatrixXd&,Eigen::VectorXi&,Eigen::VectorXd&)> nearest_neighbor = [](
    const Eigen::MatrixXd & X,
    const Eigen::MatrixXd & Y,
    Eigen::VectorXi & I,
    Eigen::VectorXd & sqrD)
  {
    const int k = X.rows();
    const int n = Y.rows();
    if(Y.cols() == 2 || Y.cols() == 3)
    {
      // O(k log k + n log k)
      Eigen::MatrixXd _;
      Eigen::VectorXi Ele = igl::LinSpaced<Eigen::VectorXi>(k,0,k-1);
      igl::point_mesh_squared_distance(Y,X,Ele,sqrD,I,_);
    }else
    {
      // O(nk)
      sqrD.setConstant(n,1,std::numeric_limits<double>::infinity());
      I.setConstant(n,1,-1);
      igl::parallel_for(n,[&](const int j)
      {
        for(int i = 0;i<k;i++)
        {

          const double sqrdji = (X.row(i)-Y.row(j)).squaredNorm();
          if(sqrdji < sqrD(j))
          {
            sqrD(j) = sqrdji;
            I(j) = i;
          }
        }
      },10000);
    }
  };
  return kmeans(P,W,k,nearest_neighbor,C,I);
}
