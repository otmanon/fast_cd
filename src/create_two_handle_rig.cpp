#include "create_two_handle_rig.h"
#include "HandleRig.h"
#ifdef WIN32
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
void create_two_handle_rig(std::string file_path, Eigen::MatrixXd& X)
{
#ifdef WIN32
    namespace fs = std::filesystem;
#else
    namespace fs = std::experimental::filesystem;
#endif
    if (!fs::exists(fs::path(file_path).parent_path()))
    {
        fs::create_directories(fs::path(file_path).parent_path());
    }

    std::vector<Eigen::Matrix4f> P(2);
    P[0].setIdentity(); P[1].setIdentity();

    
    Eigen::MatrixXd W;
    W.resize(X.rows(), 2);
    W.setZero();

    double minX = X.col(0).minCoeff();
    double maxX = X.col(0).maxCoeff();
    double width = (maxX - minX);

    //get the 10% left most vertices
    double cutoff_min = minX + width * 0.1;
    double cutoff_max = maxX - width * 0.1;

    Eigen::RowVector3d left_mean, right_mean;
    left_mean.setZero(); right_mean.setZero();

    int n_left = 0;
    int n_right = 0;
    for (int i = 0; i < X.rows(); i++)
    {
        if (X(i, 0) < cutoff_min)
        {
            W(i, 0) = 1;
            left_mean += X.row(i);
            n_left += 1;
        }
        if (X(i, 0) > cutoff_max)
        {
            W(i, 1) = 1;
            right_mean += X.row(i);
            n_right += 1;
        }
    }
    left_mean /= n_left;
    right_mean /= n_right;
    // igl::unique(F, surfaceI);

    P[0].block(0, 3, 3, 1) = left_mean.transpose().cast<float>();
    P[1].block(0, 3, 3, 1) = right_mean.transpose().cast<float>();
    HandleRig rig = HandleRig(X, P, W);
    rig.write_rig_to_json(file_path);
}