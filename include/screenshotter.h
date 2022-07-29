#pragma once
#include <Eigen/Core>
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
#include <filesystem>
struct Screenshotter
{
    char name[8] = "";
    void grab(int step, igl::opengl::glfw::Viewer* v, std::string dir, int id)
    {
        int width = v->core().viewport(2);
        int height = v->core().viewport(3);

        Eigen::Matrix<unsigned char, -1, -1> R(width, height);
        Eigen::Matrix<unsigned char, -1, -1> G(width, height);
        Eigen::Matrix<unsigned char, -1, -1> B(width, height);//
        Eigen::Matrix<unsigned char, -1, -1> A(width, height);
        //
         // Draw the scene in the buffers
        v->core().draw_buffer(v->data_list[id], false, R, G, B, A);

        sprintf(name, "%04i.png", step);
        std::string filepath = dir + "/" + name;

        if (!std::filesystem::exists(std::filesystem::path(filepath).parent_path()))
        {
            std::filesystem::create_directories(std::filesystem::path(filepath).parent_path());
        }
        // A.setOnes();
         // Save it to a PNG
        igl::png::writePNG(R, G, B, A, filepath);
    }

    void grab(std::string name, igl::opengl::glfw::Viewer* v, std::string dir, int id)
    {
        int width = v->core().viewport(2);
        int height = v->core().viewport(3);

        Eigen::Matrix<unsigned char, -1, -1> R(width, height);
        Eigen::Matrix<unsigned char, -1, -1> G(width, height);
        Eigen::Matrix<unsigned char, -1, -1> B(width, height);//
        Eigen::Matrix<unsigned char, -1, -1> A(width, height);
        //
         // Draw the scene in the buffers
        v->core().draw_buffer(v->data_list[id], false, R, G, B, A);
        std::string filepath = dir + "/" + name + ".png";

        if (!std::filesystem::exists(std::filesystem::path(filepath).parent_path()))
        {
            std::filesystem::create_directories(std::filesystem::path(filepath).parent_path());
        }
        // A.setOnes();
         // Save it to a PNG
        igl::png::writePNG(R, G, B, A, filepath);
    }

};