#pragma once
#include "fast_cd_viewer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/bind_vertex_attrib_array.h>
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

    void launch(bool launch_rendering = true)
    { //TODO initialized only data(0)
     //   fast_cd_viewer::launch();
     //   return;
        if (igl_v->data().meshgl.is_initialized)
        {
            return;
        }
        igl_v->data().meshgl.is_initialized = true;
        std::string mesh_vertex_shader_string = v_sh;
            //R"(#version 150
            //    uniform mat4 view;
            //    uniform mat4 proj;
            //    uniform mat4 normal_matrix;
            //    in vec3 position;
            //    in vec3 normal;
            //    out vec3 position_eye;
            //    out vec3 normal_eye;
            //    in vec4 Ka;
            //    in vec4 Kd;
            //    in vec4 Ks;
            //    in vec2 texcoord;
            //    out vec2 texcoordi;
            //    out vec4 Kai;
            //    out vec4 Kdi;
            //    out vec4 Ksi;

            //    void main()
            //    {
            //        position_eye = vec3 (view * vec4 (position, 1.0));
            //        normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
            //        normal_eye = normalize(normal_eye);
            //        gl_Position = proj * vec4 (position_eye, 1.0); //proj * view * vec4(position, 1.0);"
            //        Kai = Ka;
            //        Kdi = Kd;
            //        Ksi = Ks;
            //        texcoordi = texcoord;
            //    }
            // )";

            std::string mesh_fragment_shader_string =  f_sh;
 
    std::string overlay_vertex_shader_string =
    R"(#version 150
        uniform mat4 view;
        uniform mat4 proj;
        in vec3 position;
        in vec3 color;
        out vec3 color_frag;

        void main()
        {
            gl_Position = proj * view * vec4 (position, 1.0);
            color_frag = color;
        }
        )";

    std::string overlay_fragment_shader_string =
    R"(#version 150
        in vec3 color_frag;
        out vec4 outColor;
        void main()
        {
        outColor = vec4(color_frag, 1.0);
        }
        )";

    std::string overlay_point_fragment_shader_string =
        R"(#version 150
        in vec3 color_frag;
        out vec4 outColor;
        void main()
        {
            if (length(gl_PointCoord - vec2(0.5)) > 0.5)
            discard;
            outColor = vec4(color_frag, 1.0);
        }
        )";

    std::string text_vert_shader =
        R"(#version 330
        in vec3 position;
        in float character;
        in float offset;
        uniform mat4 view;
        uniform mat4 proj;
        out int vCharacter;
        out float vOffset;
        void main()
        {
            vCharacter = int(character);
            vOffset = offset;
            gl_Position = proj * view * vec4(position, 1.0);
        }
        )";

    std::string text_geom_shader =
    R"(#version 150 core
        layout(points) in;
        layout(triangle_strip, max_vertices = 4) out;
        out vec2 gTexCoord;
        uniform mat4 view;
        uniform mat4 proj;
        uniform vec2 CellSize;
        uniform vec2 CellOffset;
        uniform vec2 RenderSize;
        uniform vec2 RenderOrigin;
        uniform float TextShiftFactor;
        in int vCharacter[1];
        in float vOffset[1];
        void main()
        {
            // Code taken from https://prideout.net/strings-inside-vertex-buffers
            // Determine the final quad's position and size:
            vec4 P = gl_in[0].gl_Position + vec4( vOffset[0]*TextShiftFactor, 0.0, 0.0, 0.0 ); // 0.04
            vec4 U = vec4(1, 0, 0, 0) * RenderSize.x; // 1.0
            vec4 V = vec4(0, 1, 0, 0) * RenderSize.y; // 1.0

            // Determine the texture coordinates:
            int letter = vCharacter[0]; // used to be the character
            letter = clamp(letter - 32, 0, 96);
            int row = letter / 16 + 1;
            int col = letter % 16;
            float S0 = CellOffset.x + CellSize.x * col;
            float T0 = CellOffset.y + 1 - CellSize.y * row;
            float S1 = S0 + CellSize.x - CellOffset.x;
            float T1 = T0 + CellSize.y;

            // Output the quad's vertices:
            gTexCoord = vec2(S0, T1); gl_Position = P - U - V; EmitVertex();
            gTexCoord = vec2(S1, T1); gl_Position = P + U - V; EmitVertex();
            gTexCoord = vec2(S0, T0); gl_Position = P - U + V; EmitVertex();
            gTexCoord = vec2(S1, T0); gl_Position = P + U + V; EmitVertex();
            EndPrimitive();
        }
    )";

    std::string text_frag_shader =
    R"(#version 330
        out vec4 outColor;
        in vec2 gTexCoord;
        uniform sampler2D font_atlas;
        uniform vec3 TextColor;
        void main()
        {
            float A = texture(font_atlas, gTexCoord).r;
            outColor = vec4(TextColor, A);
        }
        )";
    igl_v->launch_init(true, false, "fast CD App", 1920, 1080);
   // igl_v->data().meshgl.init_buffers();
    init_buffers(0);
    igl_v->data().meshgl.init_text_rendering();
    igl::opengl::create_shader_program(
        mesh_vertex_shader_string,
        mesh_fragment_shader_string,
        {},
        igl_v->data().meshgl.shader_mesh);
    igl::opengl::create_shader_program(
        overlay_vertex_shader_string,
        overlay_fragment_shader_string,
        {},
        igl_v->data().meshgl.shader_overlay_lines);
    igl::opengl::create_shader_program(
        overlay_vertex_shader_string,
        overlay_point_fragment_shader_string,
        {},
        igl_v->data().meshgl.shader_overlay_points);
    igl::opengl::create_shader_program(
        text_geom_shader,
        text_vert_shader,
        text_frag_shader,
        {},
        igl_v->data().meshgl.shader_text);

    if (launch_rendering)
        {
            igl_v->core().animation_max_fps =2000;
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
        // Mesh: Vertex Array Object & Buffer objects
        glGenVertexArrays(1, &g.vao_mesh);
        glBindVertexArray(g.vao_mesh);
        glGenBuffers(1, &g.vbo_V);
        glGenBuffers(1, &g.vbo_V_normals);
        glGenBuffers(1, &g.vbo_V_ambient);
        glGenBuffers(1, &g.vbo_V_diffuse);
        glGenBuffers(1, &g.vbo_V_specular);
        glGenBuffers(1, &g.vbo_V_uv);
        glGenBuffers(1, &g.vbo_F);
        for (int i = 0; i < num_primary_vec4; i++)
        {
            glGenBuffers(1, &vbo_p_W_list[i]);
        }
        for (int i = 0; i < num_secondary_vec4; i++)
        {
            glGenBuffers(1, &vbo_s_W_list[i]);
        }
        /*   glGenBuffers(1, &vbo_pW2);
        glGenBuffers(1, &vbo_pW3);
        glGenBuffers(1, &vbo_pW4);*/
        // glGenBuffers(1, &vbo_sW);             //generate the buffers
        glGenTextures(1, &g.vbo_tex);
        glGenTextures(1, &g.font_atlas);

        // Line overlay
        glGenVertexArrays(1, &g.vao_overlay_lines);
        glBindVertexArray(g.vao_overlay_lines);
        glGenBuffers(1, &g.vbo_lines_F);
        glGenBuffers(1, &g.vbo_lines_V);
        glGenBuffers(1, &g.vbo_lines_V_colors);

        // Point overlay
        glGenVertexArrays(1, &g.vao_overlay_points);
        glBindVertexArray(g.vao_overlay_points);
        glGenBuffers(1, &g.vbo_points_F);
        glGenBuffers(1, &g.vbo_points_V);
        glGenBuffers(1, &g.vbo_points_V_colors);

        // Text Labels
        g.vertex_labels.init_buffers();
        g.face_labels.init_buffers();
        g.custom_labels.init_buffers();

        g.dirty = igl::opengl::MeshGL::DIRTY_ALL;
            
    }
    //destroy all the vbos
    void free_buffers(int id)
    {
        igl::opengl::MeshGL& g = igl_v->data_list[id].meshgl;
        if (g.is_initialized)
        {
            glDeleteVertexArrays(1, &g.vao_mesh);
            glDeleteVertexArrays(1, &g.vao_overlay_lines);
            glDeleteVertexArrays(1, &g.vao_overlay_points);

            glDeleteBuffers(1, &g.vbo_V);
            glDeleteBuffers(1, &g.vbo_V_normals);
            glDeleteBuffers(1, &g.vbo_V_ambient);
            glDeleteBuffers(1, &g.vbo_V_diffuse);
            glDeleteBuffers(1, &g.vbo_V_specular);
            glDeleteBuffers(1, &g.vbo_V_uv);
            for (int i = 0; i < num_primary_vec4; i++)
            {
                glDeleteBuffers(1, &vbo_p_W_list[i]);
            }
            for (int i = 0; i < num_secondary_vec4; i++)
            {
                glDeleteBuffers(1, &vbo_s_W_list[i]);
            }
            //  glDeleteBuffers(1, &vbo_pW1);      
        /*     glDeleteBuffers(1, &vbo_pW2);
            glDeleteBuffers(1, &vbo_pW3);
            glDeleteBuffers(1, &vbo_pW4);*/
            // glDeleteBuffers(1, &vbo_sW);
            glDeleteBuffers(1, &g.vbo_F);
            glDeleteBuffers(1, &g.vbo_lines_F);
            glDeleteBuffers(1, &g.vbo_lines_V);
            glDeleteBuffers(1, &g.vbo_lines_V_colors);
            glDeleteBuffers(1, &g.vbo_points_F);
            glDeleteBuffers(1, &g.vbo_points_V);
            glDeleteBuffers(1, &g.vbo_points_V_colors);

            // Text Labels
            g.vertex_labels.free_buffers();
            g.face_labels.free_buffers();
            g.custom_labels.free_buffers();

            glDeleteTextures(1, &g.vbo_tex);
            glDeleteTextures(1, &g.font_atlas);
        }
    }
    
    //stores the data of the primary weights and sets dirty flag to true. DOES NOT BIND THESE WEIHTS TO THE VBO
    void set_primary_weights(const MatrixXd& pW, int fid)
    {
        //ensure pW has enough colums
        RowMatrixXf p_W = pW.cast<float>();
        //s_W = sW.cast<float>();
        if (p_W.cols() > max_num_primary_bones)
        {
            printf("num primary weights needs to be less than or equal to %i\n", max_num_primary_bones);
            return;
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
        //ensure s_W has enough columns
        RowMatrixXf s_W = sW.cast<float>();
        if (s_W.cols() > max_num_secondary_bones)
        {
            printf("num secondary bones needs to be less than or equal to %i \n", max_num_secondary_bones);
            return;
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
        igl::opengl::MeshGL& g = igl_v->data_list[fid].meshgl;
        glBindVertexArray(g.vao_mesh);  //are these necessary?
        glUseProgram(g.shader_mesh);

        if (dirty_primary_weights)
        {
            for (int i = 0; i < num_primary_vec4; i++)
            {
                string name = "primary_weights_" + std::to_string(i + 1);
                igl::opengl::bind_vertex_attrib_array(g.shader_mesh, name.c_str(), vbo_p_W_list[i], p_W_list[i], true);
            }
            dirty_primary_weights = false;
        }
        if (dirty_secondary_weights)
        {
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
            dirty_primary_bones = false;
        }
        if (dirty_secondary_bones)
        {
            GLint zw = glGetUniformLocation(g.shader_mesh, "secondary_bones");
            glUniformMatrix4x3fv(zw, max_num_secondary_bones, GL_FALSE, BZ.data());
            dirty_secondary_bones = false;
        }
    }
};