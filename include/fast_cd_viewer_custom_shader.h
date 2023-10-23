#pragma once
#include "fast_cd_viewer.h"
#include "per_corner_helper.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/bind_vertex_attrib_array.h>
#include <igl/opengl/destroy_shader_program.h>


typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXf;

/*
Container for each mesh containing the fastCD info needed to render it
*/
struct fast_cd_gl
{


    int max_num_primary_bones;//number of primary bones
    int num_primary_vec4;     //number of primary vec4 we will load our weights in. (max_bones/4)

    int max_num_secondary_bones; //number of secondary bones
    int num_secondary_vec4;      //number of secondary vec4 we will load our weights in. 


 
    std::vector<RowMatrixXf> p_W_list;
    std::vector<GLuint> vbo_p_W_list;

    std::vector<RowMatrixXf> s_W_list;
    std::vector<GLuint> vbo_s_W_list;

    RowMatrixXf BP; //primary bone matrices
    RowMatrixXf BZ; //secondary bone matrices

    bool dirty_primary_weights;
    bool dirty_secondary_weights;

    bool dirty_primary_bones;
    bool dirty_secondary_bones;

    int id; //if to find this in the data_list

    fast_cd_gl(int id, int max_num_primary_bones=16, int max_num_secondary_bones=16)
    {
        this->id = id;
        this->max_num_primary_bones = max_num_primary_bones;
        this->max_num_secondary_bones = max_num_primary_bones;

        num_primary_vec4 = max_num_primary_bones / 4;
        num_secondary_vec4 = max_num_secondary_bones / 4;

        vbo_p_W_list.resize(num_primary_vec4);
        vbo_s_W_list.resize(num_secondary_vec4);

        dirty_primary_weights = true;
        dirty_secondary_weights = true;

        dirty_primary_bones = true;
        dirty_secondary_bones = true;
    }

    void init_buffers(igl::opengl::MeshGL& g)
    {
        GLint h;
        glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &h);
        glBindVertexArray(g.vao_mesh);

        for (int i = 0; i < num_primary_vec4; i++)
        {
            glGenBuffers(1, &vbo_p_W_list[i]);
        }
        for (int i = 0; i < num_secondary_vec4; i++)
        {
            glGenBuffers(1, &vbo_s_W_list[i]);//
        }

        g.dirty = igl::opengl::MeshGL::DIRTY_ALL;

        dirty_primary_weights = true;
        dirty_secondary_weights = true;
        glBindVertexArray(h);
    }

    void free_buffers(igl::opengl::MeshGL& g)
    {
        GLint h;
        glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &h);
        glBindVertexArray(g.vao_mesh);
        if (g.is_initialized)
        {

            for (int i = 0; i < num_primary_vec4; i++)
            {
                glDeleteBuffers(1, &vbo_p_W_list[i]);
            }
            for (int i = 0; i < num_secondary_vec4; i++)
            {
                glDeleteBuffers(1, &vbo_s_W_list[i]);
            }
        }
        glBindVertexArray(h);
    }
    void set_primary_weights(RowMatrixXf& p_W)
    {
        if (p_W.cols() > max_num_primary_bones)
        {
            printf("num primary weights needs to be less than or equal to %i \n", max_num_primary_bones);
            // return;
        }
        //pad with 0 if not enough columns
        if (p_W.cols() < max_num_primary_bones)
        {
            int num_cols = p_W.cols();
            p_W.conservativeResize(p_W.rows(), max_num_primary_bones);
            p_W.rightCols(max_num_primary_bones - num_cols).setZero();
        }
        for (int i = 0; i < num_primary_vec4; i++)
        {
            RowMatrixXf W = p_W.block(0, 4 * i, p_W.rows(), 4);
            p_W_list.push_back(W);
        }
        dirty_primary_weights = true;

    }

    void set_secondary_weights(RowMatrixXf& s_W)
    {
        if (s_W.cols() > max_num_secondary_bones)
        {
            printf("num secondary bones needs to be less than or equal to %i \n", max_num_secondary_bones);
            //    return;
        }
        //pad with 0 if not enough columns
        if (s_W.cols() < max_num_secondary_bones)
        {
            int num_cols = s_W.cols();
            s_W.conservativeResize(s_W.rows(), max_num_secondary_bones);
            s_W.rightCols(max_num_secondary_bones - num_cols).setZero();
        }

        for (int i = 0; i < num_secondary_vec4; i++)
        {
            RowMatrixXf W = s_W.block(0, 4 * i, s_W.rows(), 4);
            s_W_list.push_back(W);
        }
        dirty_secondary_weights = true;

    }

    void set_primary_bone_transforms(VectorXd& p)
    {
        VectorXf pf = p.cast<float>();
        MatrixXf P = Map<MatrixXf>(pf.data(), pf.rows() / 3, 3);
        if (P.rows() < 4 * max_num_primary_bones)
        {
            //pad with zeros
            int num_rows = P.rows();
            int max_num_rows = max_num_primary_bones * 4;
            P.conservativeResize(max_num_rows, 3);
            P.bottomRows(max_num_rows - num_rows).setZero();
        }
        else if (P.rows() > 4 * max_num_primary_bones)
        {
            P = P.topRows(4 * max_num_primary_bones);
        }
        BP = P;
        dirty_primary_bones = true;
    }

    void set_secondary_bone_transforms(VectorXd& z)
    {
        VectorXf zf = z.cast<float>();
        MatrixXf Z = Map<MatrixXf>(zf.data(), zf.rows() / 3, 3);
        if (Z.rows() < 4 * max_num_secondary_bones)
        {
            //pad with zeros
            int num_rows = Z.rows();
            int max_num_rows = max_num_secondary_bones * 4;
            Z.conservativeResize(max_num_rows, 3);
            Z.bottomRows(max_num_rows - num_rows).setZero();
        }
        else if (Z.rows() > 4 * max_num_secondary_bones)
        {
            Z = Z.topRows(4 * max_num_secondary_bones);
        }
        BZ = Z;
        dirty_secondary_bones = true;
    }
    void updateGL(igl::opengl::MeshGL& g)
    {
        GLint h;
        glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &h);
        glBindVertexArray(g.vao_mesh);
        //  glBindVertexArray(g.vao_mesh);  //are these necessary?//
         //glVertex
         glUseProgram(g.shader_mesh);
        //  g.v
        if (dirty_primary_weights)
        {
//            printf("Dirty primary weights \n");
            // glBindVertexBuffer(vertexBindingPoint, mesh.vbo, mesh.vboOffset, sizeof(Vertex));
            for (int i = 0; i < num_primary_vec4; i++)
            {
                string name = "primary_weights_" + std::to_string(i + 1);
                igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_p_W_list[i], p_W_list[i], true);
            }
            dirty_primary_weights = false;
        }
        if (dirty_secondary_weights)
        {
//            printf("Dirty primary weights \n");
            for (int i = 0; i < num_secondary_vec4; i++)//
            {
                string name = "secondary_weights_" + std::to_string(i + 1);
                igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_s_W_list[i], s_W_list[i], true);
            }
            dirty_secondary_weights = false;
        }
        if (dirty_primary_bones)
        {
            GLint pw = glGetUniformLocation(g.shader_mesh, "primary_bones");
            glUniformMatrix4x3fv(pw, max_num_primary_bones, GL_FALSE, BP.data());
            dirty_primary_bones =  false;
        }
        if (dirty_secondary_bones)
        {
            GLint zw = glGetUniformLocation(g.shader_mesh, "secondary_bones");
            glUniformMatrix4x3fv(zw, max_num_secondary_bones, GL_FALSE, BZ.data());
            dirty_secondary_bones = false;
        }
        glBindVertexArray(h);
    }
};

using namespace std;
struct fast_cd_viewer_custom_shader : public fast_cd_viewer
{

    
    string v_sh; //vertex shader string
    string f_sh; //fragment shader string

    int max_num_primary_bones;//number of primary bones
    int num_primary_vec4;     //number of primary vec4 we will load our weights in. (max_bones/4)

    int max_num_secondary_bones; //number of secondary bones
    int num_secondary_vec4;      //number of secondary vec4 we will load our weights in. 

    std::vector<fast_cd_gl> fcd_gl;


    fast_cd_viewer_custom_shader() {};
	fast_cd_viewer_custom_shader(string& vertex_shader, string& fragment_shader, int max_b_p=16, int max_b_s = 16) :fast_cd_viewer()
	{
        namespace fs = std::filesystem;
        if (!fs::exists(fs::path(vertex_shader)))
        {
            printf("Could not find vertex shader!\n");
            return;
        }
        if (!fs::exists(fs::path(fragment_shader)))
        {
            printf("Coult not find fragment shader\n"); //
            return;
        }
        ifstream f_v(vertex_shader); //taking file as inputstream

        if (f_v) {
            ostringstream ss;
            ss << f_v.rdbuf(); // reading data
            v_sh = ss.str();
        }

        ifstream f_f(fragment_shader); //taking file as inputstream
    
        if (f_f) {
            ostringstream ss;
            ss << f_f.rdbuf(); // reading data
            f_sh = ss.str();
        }
       
       max_num_secondary_bones = max_b_s;
       max_num_primary_bones = max_b_p;
       if (max_num_primary_bones % 4 != 0)
           printf("Maximum number of primary bones must be a multiple of 4 !\n");
       if (max_num_secondary_bones % 4 != 0)
           printf("Maximum number of secondary bones must be a multiple of 4 !\n");
       num_primary_vec4 = max_num_primary_bones / 4;
       num_secondary_vec4 = max_num_secondary_bones / 4;
       fast_cd_gl f = fast_cd_gl(0, max_num_primary_bones, max_num_primary_bones);

       fcd_gl.push_back(f);
    }


    void init_all_shaders()
    {
        for (int i = 0; i < igl_v->data_list.size(); i++)
        {
            init_shaders(i);
        }
    }

    void init_shaders(int id)
    {
        igl_v->data_list[id].meshgl.free();
        igl_v->data_list[id].meshgl.init();

        init_buffers(id);


        igl::opengl::MeshGL& g = igl_v->data_list[id].meshgl;
        igl::opengl::destroy_shader_program(g.shader_mesh);

        std::string mesh_vertex_shader_string = v_sh;

        std::string mesh_fragment_shader_string = f_sh;

        igl::opengl::create_shader_program(
            mesh_vertex_shader_string,
            mesh_fragment_shader_string,
            {},
            g.shader_mesh);

    }
    void launch(int max_fps = 6000, bool launch_rendering = true)
    { 

    igl_v->launch_init(true, false, "fast CD App", 1920, 1080);

    init_all_shaders();
    if (launch_rendering)
        {
        igl_v->core().animation_max_fps = max_fps;// max_fps;
            igl_v->core().is_animating = true;
            igl_v->launch_rendering(true);
            igl_v->launch_shut();
            //fast_cd_viewer::launch();
           // igl_v->launch(true, false, "fast cd app", 1920, 1080);
        }
    }

    //vreate all the vbo indicies
    void init_buffers(int id)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[id].meshgl;
        assert(id < igl_v->data_list.size() && "id is not associated with an igl data mesh!");
        fcd_gl[id].init_buffers(g);
        g.dirty = igl::opengl::MeshGL::DIRTY_ALL;
    }
    //destroy all the vbos
    void free_buffers(int id)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[id].meshgl;
        fcd_gl[id].free_buffers(g);
    }
    
    //stores the data of the primary weights and sets dirty flag to true. DOES NOT BIND THESE WEIHTS TO THE VBO
    void set_primary_weights(const MatrixXd& pW, int fid)
    {
 
        RowMatrixXf p_W;
        if (igl_v->data_list[fid].face_based)
        {
            per_corner_helper(pW, igl_v->data_list[fid].F, p_W);
        }
        else
        {
            //ensure pW has enough colums
             p_W = pW.cast<float>();
        }

        fcd_gl[fid].set_primary_weights(p_W);

    }

    //stores the data of the secondary weights and sets dirty flag to true. DOES NOT BIND THESE WEIHTS TO THE VBO
    void set_secondary_weights(const MatrixXd& sW, int fid)
    { 
   
        RowMatrixXf s_W;
        if (igl_v->data_list[fid].face_based)
        {
            per_corner_helper(sW, igl_v->data_list[fid].F, s_W);
        }
        else
        {
            //ensure pW has enough colums
            s_W = sW.cast<float>();
        }
        
        fcd_gl[fid].set_secondary_weights(s_W);
  
    }

    void set_weights(const MatrixXd& pW, const MatrixXd& sW, int fid)
    {
        
        set_primary_weights(pW, fid);
        set_secondary_weights(sW, fid);
    }

    void set_primary_bone_transforms(VectorXd& p, int id)
    {
    
        fcd_gl[id].set_primary_bone_transforms(p);

    }

    void set_secondary_bone_transforms(VectorXd& z, int id)
    {

        fcd_gl[id].set_secondary_bone_transforms(z);

    }

    void set_bone_transforms(VectorXd& p, VectorXd& z, int id)
    {
   
        set_primary_bone_transforms(p, id);
        set_secondary_bone_transforms(z, id);
     

    }

   
    void add_mesh(int& id)
    {
        igl_v->append_mesh();
        id = igl_v->data_list.size() - 1;
        igl_v->data_list[id].clear();
   //     init_all_shaders();
        fast_cd_gl f = fast_cd_gl(id, max_num_primary_bones, max_num_secondary_bones);
        fcd_gl.push_back(f);
    }

    void add_mesh(const MatrixXd& V, const MatrixXi& F, int& id)
    {
        igl_v->append_mesh();
        id = igl_v->data_list.size() - 1;
        igl_v->data_list[id].clear();

    //    init_all_shaders();
        fast_cd_gl f = fast_cd_gl(id, max_num_primary_bones, max_num_secondary_bones);
        igl_v->data_list[id].set_mesh(V, F);

        fcd_gl.push_back(f);
    }


    //call this to update any of the uniforms/vertex attributes that have changes by checking the dirty flags
    void updateGL(int fid)
    {
        igl_v->data_list[fid].face_based = false; // temporary fix, make sure face based is false
        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
        fcd_gl[fid].updateGL(g);

    }
};