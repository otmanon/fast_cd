#include "line_search.h"
#include <iostream>

double line_search_pgd(
    const std::function<double(const Eigen::VectorXd&)>& f,
    const std::function<void(Eigen::VectorXd&)>& proj_z,
    const Eigen::VectorXd& z,
    const Eigen::VectorXd& dz,
    const double max_step)
{
    // Huh. This is not "by the book"
    // https://en.wikipedia.org/wiki/Backtracking_line_search 
    // Should this really include the projection? If not, then the step might
    // decrease the energy but then the projection makes things worse... If the
    // projection is included then we're not really searching along a line.
    const double m = -1e-10;
    const double tao = 0.5;
    const double c = 0.5;
    double t = -c * m;
    const int max_iters = 100;
    const double f0 = f(z);
    double alpha = max_step;
   
    for (int j = 0; j < max_iters; j++)
    {
        Eigen::VectorXd zj = z + alpha * dz;
        proj_z(zj);
        const double fj = f(zj);
        if ((f0 - fj) >= alpha * t)
        {
            break;
        }
        alpha *= tao;
    }
    return alpha;

}


double line_search(const std::function<double(const Eigen::VectorXd&)>& f,
    const Eigen::VectorXd& z,
    const Eigen::VectorXd& dz, const double max_step)
{
    const double m = -1e-10;
    const double tao = 0.5;
    const double c = 0.5;
    double t = -c * m;
    const int max_iters = 20;
    const double f0 = f(z);
    double alpha = max_step;

    for (int j = 0; j < max_iters; j++)
    {
        Eigen::VectorXd zj = z + alpha * dz;
        const double fj = f(zj);
        if ((f0 - fj) >= alpha * t)
        {
            break;
        }
        alpha *= tao;
    }
    return alpha;
}