#pragma once
#include "fast_cd_viewer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/bind_vertex_attrib_array.h>
#include <igl/opengl/destroy_shader_program.h>

//if face based
  // Input:
//   X  #V by dim quantity
// Output:
//   X_vbo  #F*3 by dim scattering per corner
void per_corner(const Eigen::MatrixXd& X, const Eigen::MatrixXi& F,
    igl::opengl::MeshGL::RowMatrixXf& X_vbo)
{
    X_vbo.resize(F.rows() * 3, X.cols());
    for (unsigned i = 0; i < F.rows(); ++i)
        for (unsigned j = 0; j < 3; ++j)
            X_vbo.row(i * 3 + j) = X.row(F(i, j)).cast<float>();
   
}


using namespace std;
struct fast_cd_viewer_custom_shader : public fast_cd_viewer
{

    string v_sh; //vertex shader string
    string f_sh; //fragment shader string

    int max_num_primary_bones;//number of primary bones
    int num_primary_vec4;     //number of primary vec4 we will load our weights in. (max_bones/4)

    int max_num_secondary_bones; //number of secondary bones
    int num_secondary_vec4;      //number of secondary vec4 we will load our weights in. 


    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXf;
    std::vector<RowMatrixXf> p_W_list;
    std::vector<GLuint> vbo_p_W_list;

    std::vector<RowMatrixXf> s_W_list;
    std::vector<GLuint> vbo_s_W_list;

    RowMatrixXf BP; //primary bone matrices
    RowMatrixXf BZ; //secondary bone matrices
    //secondary weights and primary weights will always be dirty at the same time: at the star tof the sim. For some reason I can't initialize them until the first draw call
  
    bool dirty_primary_weights;
    bool dirty_secondary_weights;

    bool dirty_primary_bones;
    bool dirty_secondary_bones;

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
       printf("Assuming maximum number of bones is %i, please make sure vertex \
            buffer has this value set in the primary_bones[n] uniform, where n==max_num_bones\n", max_num_primary_bones);
       num_primary_vec4 = max_num_primary_bones / 4;
       num_secondary_vec4 = max_num_secondary_bones / 4;
        //initialize list 
       vbo_p_W_list.resize(num_primary_vec4);
       vbo_s_W_list.resize(num_secondary_vec4);

       dirty_primary_weights = false;
       dirty_secondary_weights = false;

       dirty_primary_bones = false;
       dirty_secondary_bones = false;
    }

    void launch(int max_fps = 60, bool launch_rendering = true)
    { 

    igl_v->launch_init(true, false, "fast CD App", 1920, 1080);

    igl_v->data_list[0].meshgl.init();
 
    init_buffers(0);

   
    igl::opengl::MeshGL& g = igl_v->data_list[0].meshgl;
    igl::opengl::destroy_shader_program(g.shader_mesh);

    std::string mesh_vertex_shader_string = v_sh;

    std::string mesh_fragment_shader_string =  f_sh;

    igl::opengl::create_shader_program(
        mesh_vertex_shader_string,
        mesh_fragment_shader_string,
        {},
        g.shader_mesh);

    if (launch_rendering)
        {
            igl_v->core().animation_max_fps = max_fps;
            igl_v->launch_rendering(true);
            igl_v->launch_shut();
            //fast_cd_viewer::launch();
           // igl_v->launch(true, false, "fast cd app", 1920, 1080);
        }
    }

    //vreate all the vbo indicies
    void init_buffers(int id)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[0].meshgl;
      
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
            
    }
    //destroy all the vbos
    void free_buffers(int id)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[id].meshgl;
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
    }
    


    //stores the data of the primary weights and sets dirty flag to true. DOES NOT BIND THESE WEIHTS TO THE VBO
    void set_primary_weights(const MatrixXd& pW, int fid)
    {

        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
 
        RowMatrixXf p_W;
        if (igl_v->data_list[fid].face_based)
        {
            per_corner(pW, igl_v->data_list[fid].F, p_W);
        }
        else
        {
            //ensure pW has enough colums
             p_W = pW.cast<float>();
        }

        if (p_W.cols() > max_num_primary_bones)
        {
            printf("num primary weights needs to be less than or equal to %i\n", max_num_primary_bones);
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

    //stores the data of the secondary weights and sets dirty flag to true. DOES NOT BIND THESE WEIHTS TO THE VBO
    void set_secondary_weights(const MatrixXd& sW, int fid)
    { 

        RowMatrixXf s_W;
        if (igl_v->data_list[fid].face_based)
        {
            per_corner(sW, igl_v->data_list[fid].F, s_W);
        }
        else
        {
            //ensure pW has enough colums
            s_W = sW.cast<float>();
        }
        
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

    void set_weights(const MatrixXd& pW, const MatrixXd& sW, int fid)
    {
        set_primary_weights(pW, fid);
        set_secondary_weights(sW, fid);
    }

    void set_primary_bone_transforms(VectorXd& p, int id)
    {
        VectorXf pf = p.cast<float>();
        MatrixXf P = Map<MatrixXf>(pf.data(), p.rows() / 3, 3);
        if (P.rows() < 4 * max_num_primary_bones)
        {
            //pad with zeros
            int num_rows = P.rows();
            int max_num_rows = max_num_primary_bones * 4;
            P.conservativeResize(max_num_rows, 3);
            P.bottomRows(max_num_rows - num_rows).setZero();
            pf = Map<VectorXf>(P.data(), P.rows() * 3, 1);
        }
        else if (P.rows() > 4 * max_num_primary_bones)
        {
            P = P.topRows(4 * max_num_primary_bones);
            pf = Map<VectorXf>(P.data(), P.rows() * 3, 1);
        }
        BP = P;
        dirty_primary_bones = true;
    }

    void set_secondary_bone_transforms(VectorXd& z, int id)
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
            zf = Map<VectorXf>(Z.data(), Z.rows() * 3, 1);
        }
        else if (Z.rows() > 4 * max_num_secondary_bones)
        {
            Z = Z.topRows(4 * max_num_secondary_bones);
            zf = Map<VectorXf>(Z.data(), Z.rows() * 3, 1);
        }
        BZ = Z;
        dirty_secondary_bones = true;
    }

    void set_bone_transforms(VectorXd& p, VectorXd& z, int id)
    {
        set_primary_bone_transforms(p, id);
        set_secondary_bone_transforms(z, id);

    }

    /// DEPRACATING THESE BIND FUNCTIONS 
    void bind_primary_weights(int fid)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
        glBindVertexArray(g.vao_mesh);
        glUseProgram(g.shader_mesh);
        //just set this to true, will only happen once at the start of the sim anyways
        bool dirty_weights = true;
        for (int i = 0; i < num_primary_vec4; i++)
        {
            string name = "primary_weights_" + std::to_string(i + 1);
            igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_p_W_list[i], p_W_list[i], dirty_weights);
        }
    }
    
    void bind_secondary_weights(int fid) //
    {
        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
        glBindVertexArray(g.vao_mesh);
        glUseProgram(g.shader_mesh);
        bool dirty_weights = true;
        for (int i = 0; i < num_secondary_vec4; i++)//
        {
            string name = "secondary_weights_" + std::to_string(i + 1);
            igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_s_W_list[i], s_W_list[i], dirty_weights);
        }
    }

    void bind_weights( int fid)
    {
        bind_primary_weights(fid);
        bind_secondary_weights(fid);
    }
    /*
    Binds the primary and secondary weights to our vertex buffer object.
    */
    void bind_weights(MatrixXd& pW, MatrixXd& sW, int fid)
    {
        set_primary_weights(pW, fid);
        set_secondary_weights(sW, fid);
        bind_primary_weights(fid);
        bind_secondary_weights(fid);
    }

    void bind_bone_transforms(VectorXd& p, VectorXd& z, int id)
    {
        set_bone_transforms(p, z, id);

        updateGL(id);
    }


    //call this to update any of the uniforms/vertex attributes that have changes by checking the dirty flags
    void updateGL(int fid)
    {
        igl_v->data_list[fid].face_based = false; // temporary fix, make sure face based is false
        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
      //  glBindVertexArray(g.vao_mesh);  //are these necessary?//
       //glVertex
   //    glUseProgram(g.shader_mesh);
      //  g.v
        if (dirty_primary_weights)
        {
            GLint h;
            glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &h);
           glBindVertexArray(g.vao_mesh);
           // glBindVertexBuffer(vertexBindingPoint, mesh.vbo, mesh.vboOffset, sizeof(Vertex));
            for (int i = 0; i < num_primary_vec4; i++)
            {
                string name = "primary_weights_" + std::to_string(i + 1);
                igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_p_W_list[i], p_W_list[i], true);
            }
            dirty_primary_weights = false;

         //   glBindVertexArray(h);
        }
        if (dirty_secondary_weights)   
        {
            GLint h;
           glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &h);
            glBindVertexArray(g.vao_mesh);
        //    glBindVertexArray(g.vao_mesh);
            for (int i = 0; i < num_secondary_vec4; i++)//
            {
                string name = "secondary_weights_" + std::to_string(i + 1);
                igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_s_W_list[i], s_W_list[i], true);
            }
            dirty_secondary_weights = false;
            
          //  glBindVertexArray(h);
        }
        if (dirty_primary_bones)
        {
         //   glBindVertexArray(g.vao_mesh);
            GLint pw = glGetUniformLocation(g.shader_mesh, "primary_bones");
            glUniformMatrix4x3fv(pw, max_num_primary_bones, GL_FALSE, BP.data());
            dirty_primary_bones = false;
        }
        if (dirty_secondary_bones)
        {
          //  glBindVertexArray(g.vao_mesh);
            GLint zw = glGetUniformLocation(g.shader_mesh, "secondary_bones");
            glUniformMatrix4x3fv(zw, max_num_secondary_bones, GL_FALSE, BZ.data());
            dirty_secondary_bones = false;
        }
       // glBindVertexArray(0);

    }
};