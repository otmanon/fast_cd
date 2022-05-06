#include "InteractiveCDHook.h"

#include "igl/slice.h"
#include <igl/png/readPNG.h>
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <string>
#include "rainbow_cmap.h"

#include <iostream>

#include <filesystem>

#include "create_two_handle_rig.h"

#include <igl/LinSpaced.h>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/destroy_shader_program.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/readPLY.h>
#include <igl/readMSH.h>
#include <igl/readDMAT.h>
#include <igl/colon.h>
#include <igl/slice_into.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <igl/boundary_facets.h>
//Render UI related
bool InteractiveCDHook::render(igl::opengl::glfw::Viewer& viewer)
{
  
    if (as.animation_mode == ANIMATION_MODE::INTERACTIVE_ANIMATION)
    {
        if (as.constraint_type == PINNING)
            render_pinning(viewer);
        else
            render_cd(viewer);
   
        as.rig_controller->render(viewer);
    }
    else
    {
        this->fcd_viewer.configure_matcap(as.matcap_file, V0, F);
        viewer.data_list[0].clear();
        viewer.data_list[0].invert_normals = true;
        viewer.data_list[0].double_sided = true;
        Eigen::MatrixXd U = Eigen::Map<Eigen::MatrixXd>(u_curr.data(), u_curr.rows(), u_curr.cols());
        Eigen::MatrixXd P = U + V0;
        viewer.data_list[0].set_mesh(P, F);
    }
    return false;
     
}

void InteractiveCDHook::render_cd(igl::opengl::glfw::Viewer& viewer)
{
    if (as.do_reduction)
    {
        Eigen::VectorXf z = z_curr.cast<float>();
        Eigen::VectorXf p = p_curr.cast<float>();

        if (as.proj_gpu)
        {  //send reduced parameters to GPU. Reconstruct full space there.
            render_reduced_gpu_proj(viewer, z, p, V.rows(), v_state.coarse_vis_id);
        }
        else
        {
            //construct full space in CPU, then send that to GPU
            if (v_state.vis_mode == TEXTURES)
            {
                render_reduced_cpu_proj(viewer, z, p, cd_WB, cd_WJ, v_state.fine_vis_id);
                if (v_state.show_cage)
                {
                    render_reduced_cpu_proj(viewer, z, p, cd_B_ext, cd_J_ext, v_state.coarse_vis_id);
                }
            }
            else
            {
                render_reduced_cpu_proj(viewer, z, p, cd_B_ext, cd_J_ext, v_state.coarse_vis_id);
            }
        }
    }
    else
    {   // no need to worry about proj_gpu. just use set_vertices
        Eigen::MatrixXd P_ext;
        if (v_state.vis_mode == TEXTURES)
        {

            Eigen::VectorXd r = cd_WJ * p_next;
            Eigen::MatrixXd P = Eigen::Map<Eigen::MatrixXd>(r.data(), r.rows() / 3, 3);
            if (v_state.vis_cd)
            {
                Eigen::VectorXd uc = W_low_to_high* uc_curr;
                P += Eigen::Map<Eigen::MatrixXd>(uc.data(), uc.rows()/3, 3);
            }
            igl::slice(P, ext_ind, 1, P_ext);
            render_full(viewer, P, v_state.fine_vis_id);
            if (v_state.show_cage)
            {
                render_full(viewer, P_ext, v_state.coarse_vis_id);
            }
        }
        else
        {
            Eigen::VectorXd uc = uc_curr;
            Eigen::VectorXd r = rig->J * as.rig_controller->p_rel;
            Eigen::MatrixXd P = Eigen::Map<Eigen::MatrixXd>(r.data(), r.rows() / 3, 3);
            if (v_state.vis_cd)
                P += Eigen::Map<Eigen::MatrixXd>(uc.data(), uc.rows() / 3, 3);
            igl::slice(P, ext_ind, 1, P_ext);
            render_full(viewer, P_ext, v_state.coarse_vis_id);
        }
    }
}

void InteractiveCDHook::render_pinning(igl::opengl::glfw::Viewer& viewer)
{
    if (as.do_reduction)
    {
        Eigen::VectorXf z = z_curr.cast<float>();

        if (as.proj_gpu)
        {  //send reduced parameters to GPU. Reconstruct full space there.
            render_reduced_gpu_proj_pin(viewer, z, V.rows(), v_state.coarse_vis_id);
        }
        else
        {
            //construct full space in CPU, then send that to GPU
            if (v_state.vis_mode == TEXTURES)
            {
                render_reduced_cpu_proj_pin(viewer, z,  pinned_WB, V_high_res0, v_state.fine_vis_id);
                if (v_state.show_cage)
                {
                    render_reduced_cpu_proj_pin(viewer, z,  pinned_B_ext, V0_ext, v_state.coarse_vis_id);
                }
            }
            else
            {
                render_reduced_cpu_proj_pin(viewer, z, pinned_B_ext,  V0_ext, v_state.coarse_vis_id);
            }
        }
    }
    else
    {   // no need to worry about proj_gpu. just use set_vertices

        Eigen::MatrixXd P_ext;
        if (v_state.vis_mode == TEXTURES)
        {
            Eigen::MatrixXd P = V_high_res0;
            if (v_state.vis_cd)
            {
                Eigen::VectorXd u = W_low_to_high * u_curr;
                P += Eigen::Map<Eigen::MatrixXd>(u.data(), u.rows() / 3, 3);
            }
            render_full(viewer, P, v_state.fine_vis_id);
            if (v_state.show_cage)
            {
                igl::slice(P, ext_ind, 1, P_ext);
                render_full(viewer, P_ext, v_state.coarse_vis_id);
            }
        }
        else
        {
            Eigen::MatrixXd P = V0;
            if (v_state.vis_cd)
            {
                P += Eigen::Map<Eigen::MatrixXd>(u_curr.data(), u_curr.rows() / 3, 3);
            }
            igl::slice(P, ext_ind, 1, P_ext);
            render_full(viewer, P_ext, v_state.coarse_vis_id);
        }
    }
}


void InteractiveCDHook::render_full(igl::opengl::glfw::Viewer& v, Eigen::MatrixXd& V, int cid)
{
    GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, false); //dopn't use gpu to compute displacements!!
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);

    v.data_list[cid].set_vertices(V);

}

void InteractiveCDHook::render_reduced_cpu_proj(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::VectorXf& p, Eigen::MatrixXd& B, Eigen::SparseMatrix<double>& J, int cid)
{
 
    GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, false); //dopn't use gpu to compute displacements!!
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);
    Eigen::VectorXd u =  J * p.cast<double>();  //the titular cpu proj
    if (v_state.vis_cd)
        u += B * z.cast<double>();
    Eigen::MatrixXd V = Eigen::Map<Eigen::MatrixXd>(u.data(), u.rows() / 3, 3);

    v.data_list[cid].set_vertices(V);
   
}

/*
 Sets data in the mesh by passing reduced buffers to the vertex shader.

 v - libigl viewer
 z - reduced space coefficients for our simulation
 p - rig parameters (row flattened 12x1 transformation matrices)
 n - (int) number of vertices

 */
void InteractiveCDHook::render_reduced_gpu_proj(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::VectorXf& p, int n, int cid)
{
    //TODO: assert that viewer is already configured for reduced step... otherwise need to resend texture

    const int m = z.rows();
    const int b = int(p.rows() / 12);
    const int s = ceil(sqrt(n * (m + b)));

    GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, v.data_list[cid].meshgl.vbo_tex);
    glUseProgram(prog_id);
    GLint n_loc = glGetUniformLocation(prog_id, "n");
    glUniform1i(n_loc, n);
    GLint m_loc = glGetUniformLocation(prog_id, "m");
    glUniform1i(m_loc, m);
    GLint b_loc = glGetUniformLocation(prog_id, "b");
    glUniform1i(b_loc, b);
    GLint s_loc = glGetUniformLocation(prog_id, "s");
    glUniform1i(s_loc, s);
    GLint q_loc = glGetUniformLocation(prog_id, "q");  //modal activations
    glUniform1fv(q_loc, m, z.data());
    GLint p_loc = glGetUniformLocation(prog_id, "p"); //rig parameters
    glUniform1fv(p_loc, p.rows(), p.data());
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, true); //if we call this then we must assume it's true
    // Do this now so that we can stop texture from being loaded by viewer
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);
    if (v.data_list[cid].dirty)
    {
        v.data_list[cid].updateGL(
            v.data_list[cid],
            v.data_list[cid].invert_normals,
            v.data_list[cid].meshgl
        );
        v.data_list[cid].dirty = igl::opengl::MeshGL::DIRTY_NONE;
    }

    v.data_list[cid].dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;

}

void InteractiveCDHook::render_reduced_cpu_proj_pin(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, Eigen::MatrixXd& B, Eigen::MatrixXd& X, int cid)
{
    GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, false); //dopn't use gpu to compute displacements!!
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, true);
    Eigen::VectorXd u = B * z.cast<double>();  //the titular cpu proj
    Eigen::MatrixXd V = Eigen::Map<Eigen::MatrixXd>(u.data(), u.rows() / 3, 3);
    V += X;

    v.data_list[cid].set_vertices(V);
}

void InteractiveCDHook::render_reduced_gpu_proj_pin(igl::opengl::glfw::Viewer& v, Eigen::VectorXf& z, int n, int cid)
{  
    const int m = z.rows();
    const int b = int(0);
    const int s = ceil(sqrt(n * (m + b)));

    GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, v.data_list[cid].meshgl.vbo_tex);
    GLint n_loc = glGetUniformLocation(prog_id, "n");
    glUniform1i(n_loc, n);
    GLint m_loc = glGetUniformLocation(prog_id, "m");
    glUniform1i(m_loc, m);
    GLint b_loc = glGetUniformLocation(prog_id, "b");
    glUniform1i(b_loc, b);
    GLint s_loc = glGetUniformLocation(prog_id, "s");
    glUniform1i(s_loc, s);
    GLint q_loc = glGetUniformLocation(prog_id, "q");  //modal activations
    glUniform1fv(q_loc, m, z.data());
 
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, true);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, true); //if we call this then we must assume it's true
    // Do this now so that we can stop texture from being loaded by viewer
    if (v.data_list[cid].dirty)
    {
        v.data_list[cid].updateGL(
            v.data_list[cid],
            v.data_list[cid].invert_normals,
            v.data_list[cid].meshgl
        );
        v.data_list[cid].dirty = igl::opengl::MeshGL::DIRTY_NONE;
    }

    v.data_list[cid].dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;

}

void InteractiveCDHook::init_viewer(igl::opengl::glfw::Viewer& v)
{
  
    this->viewer = &v;
    fcd_viewer = FastCDViewer(&v);

    if (viewer->data_list.size() == 1)
        this->viewer->append_mesh();
    //this->viewer->data_list[v_state.coarse_vis_id].clear();
    this->viewer->data_list[v_state.fine_vis_id].clear();
    new_v_state = v_state;

  // set_viewer_clusters();
    if (as.proj_gpu)
    {
        if (as.constraint_type == PINNING)
        {
            Eigen::MatrixXd Z;
            fcd_viewer.configure_deformation_texture(this->V0, this->F, sim.B, Z);
           // set_viewer_defo_textures(*viewer, this->V0, this->F, sim.B, Z, v_state.coarse_vis_id, v_state.coarse_vis_id);
        }
        else
            fcd_viewer.configure_deformation_texture(this->V0, this->F, sim.B, rig->W);
            //set_viewer_defo_textures(*viewer, this->V0, this->F, cd_sim.B, rig->W, v_state.coarse_vis_id, v_state.coarse_vis_id);
    }
    else
    {
        if (v_state.vis_mode == VIS_MODE::MATCAP)
        {
         //   set_viewer_matcap(*viewer, V_ext, F_ext, "../data" + as.matcap_file, v_state.coarse_vis_id, v_state.fine_vis_id);
            fcd_viewer.configure_matcap("../data/" + as.matcap_file, V_ext, F_ext);
        }
        else if (v_state.vis_mode == VIS_MODE::CLUSTERS)
        {
            Eigen::VectorXi labels_faces, labels;
            labels = as.constraint_type == PINNING ? sim.labels : cd_sim.labels;

            igl::slice(labels, FiT, labels_faces);
          //  set_viewer_clusters(*viewer, V_ext, F_ext, labels_faces, v_state.coarse_vis_id, v_state.fine_vis_id);
            fcd_viewer.configure_clusters(V_ext, F_ext, labels_faces);
        }
        else if (v_state.vis_mode == VIS_MODE::TEXTURES)
        {
            fcd_viewer.configure_color_texture(as.texture_png_path, V0_ext, F_ext, V_high_res,
                F_high_res, UV_high_res, FUV_high_res);
           // set_viewer_color_textures(*viewer, as.texture_png_path, V0_ext, F_ext, V_high_res,
           //     F_high_res, UV_high_res, FUV_high_res, v_state.coarse_vis_id, v_state.fine_vis_id);
        }
    }
}

void InteractiveCDHook::set_viewer_defo_textures(igl::opengl::glfw::Viewer& v, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W, int cid, int fid)
{
    using namespace Eigen;
    using namespace std;

    ///////////////////////////////////////////////////////////////////
    // Load and prepare data
    ///////////////////////////////////////////////////////////////////
    Eigen::Matrix< float, Eigen::Dynamic, 1> I;
    Eigen::MatrixXf U;
    {
        Eigen::MatrixXd Ud = B;
        U = Ud.cast<float>();
    }
    assert((U.rows() == V.rows() * 3) && "#U should be 3*#V");
    I = igl::LinSpaced< Eigen::Matrix< float, Eigen::Dynamic, 1> >(V.rows(), 0, V.rows() - 1);
    const int n = V.rows();
    const int m = U.cols();
    const int b = W.cols();       //number of bones
    const int s = ceil(sqrt(n * (m + b)));
    assert(s * s > n * (m + b));
    printf("verts : %d  modes : %d bones: %d im_dim %d\n", n, m, b, s);
    tex = Eigen::Matrix< float, Eigen::Dynamic, 3, Eigen::RowMajor>::Zero(s * s, 3);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            for (int c = 0; c < 3; c++)
            {
                tex(i * (m + b) + j, c) = U(i + c * n, j);
            }
        }
        for (int j = 0; j < b; j++)
        {
            for (int c = 0; c < 3; c++)
            {
                tex(i * (m + b) + j + m, c) = W(i, j); //bone is the same, place right after the modes
            }
        }
    }


    ///////////////////////////////////////////////////////////////////
    // Initialize viewer and opengl context
    ///////////////////////////////////////////////////////////////////
   // igl::opengl::glfw::Viewer v;
    
    v.data_list[cid].clear();
    v.data_list[fid].clear();
    v.data_list[cid].set_mesh(V, F);
    v.data_list[cid].set_face_based(false);
    v.data_list[cid].use_matcap = false;
    v.data_list[cid].show_lines = true;
    v.data_list[cid].show_faces = true;
    v.data_list[cid].invert_normals = false;
    v.data_list[cid].double_sided = true;
    v.data_list[cid].show_lines = false;
    v.data_list[cid].show_texture = false;
    v.data_list[fid].show_texture = false;
    v.data_list[fid].is_visible = false;
    v.data_list[cid].is_visible = true;


    v.data_list[cid].show_face_labels = false;
    v.data_list[cid].show_vertex_labels = false;

    v.data_list[fid].show_face_labels = false;
    v.data_list[fid].show_vertex_labels = false;
    ///////////////////////////////////////////////////////////////////
   // Send texture and vertex attributes to GPU... should make this its own function but maybe it's okay
   ///////////////////////////////////////////////////////////////////
    {
        GLuint prog_id = v.data_list[cid].meshgl.shader_mesh;
        glUseProgram(prog_id);
        GLuint VAO = v.data_list[cid].meshgl.vao_mesh;
        glBindVertexArray(VAO);
        GLuint IBO;
        glGenBuffers(1, &IBO);
        glBindBuffer(GL_ARRAY_BUFFER, IBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * I.size(), I.data(), GL_STATIC_DRAW);
        GLint iid = glGetAttribLocation(prog_id, "id");
        glVertexAttribPointer(
            iid, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(iid);
        //  glBindVertexArray(0);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, v.data_list[cid].meshgl.vbo_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, s, s, 0, GL_RGB, GL_FLOAT, tex.data());
    }
    //This is very important... otherwise when the new mesh is set, it triggers a dirty call , and igl overwrites the texture with junk...
    v.data_list[cid].dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;
}
void InteractiveCDHook::set_viewer_color_textures(igl::opengl::glfw::Viewer& viewer, std::string texture_filepath, Eigen::MatrixXd& V_coarse, Eigen::MatrixXi& F_coarse, Eigen::MatrixXd& V_fine, Eigen::MatrixXi& F_fine, 
    Eigen::MatrixXd& UV_fine, Eigen::MatrixXi& FUV_fine, int cid, int fid)
{

    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh
    viewer.data_list[fid].clear();
    viewer.data_list[cid].clear();

    viewer.data_list[fid].is_visible = true; 
    viewer.data_list[cid].show_lines = true;
    viewer.data_list[cid].show_faces = false;
    viewer.data_list[cid].is_visible = false;//dont wanna see or hear about you when using textures

    viewer.data_list[fid].clear();
        
    viewer.data_list[fid].invert_normals = false;
    viewer.data_list[fid].double_sided = false;
    viewer.data_list[fid].set_face_based(true);


    viewer.data_list[cid].show_face_labels = false;
    viewer.data_list[cid].show_vertex_labels = false;

    viewer.data_list[fid].show_face_labels = false;
    viewer.data_list[fid].show_vertex_labels = false;

    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string texture_path = texture_filepath;
    bool read = igl::png::readPNG(texture_path, R, G, B, A);

    viewer.data_list[cid].set_mesh(V_coarse, F_coarse);
    viewer.data_list[fid].set_mesh(V_fine, F_fine);
    if (UV_high_res.rows() > 0)  //if our texture matches the display mesh
    {
        //  viewer.data().set_normals(FN);
        viewer.data_list[fid].set_uv(UV_fine, FUV_fine);
        viewer.data_list[fid].show_lines = false;
        viewer.data_list[fid].show_texture = true;
        viewer.data_list[fid].show_faces = true;
        viewer.data_list[fid].set_texture(R, G, B);
    }else
    {
        viewer.data_list[fid].invert_normals = true;
    }
    
}


void InteractiveCDHook::set_viewer_matcap(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::string matcap_file, int cid, int fid)
{
    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh

    viewer.data_list[fid].clear();
    viewer.data_list[cid].clear();
 
    viewer.data_list[cid].set_mesh(V, F);
    viewer.data_list[fid].set_mesh(V_high_res, F_high_res); //should only send V_ext

    viewer.data_list[fid].is_visible = false; //dont wanna see or hear about you when using matcap young man... though tbh why not??
    viewer.data_list[cid].is_visible = true; //dont wanna see or hear about you when using matcap young man... though tbh why not??


    viewer.data_list[cid].point_size = 10;


    viewer.data_list[cid].show_face_labels = false;
    viewer.data_list[cid].show_vertex_labels = false;

    viewer.data_list[fid].show_face_labels = false;
    viewer.data_list[fid].show_vertex_labels = false;

    viewer.data_list[cid].invert_normals = false;
    viewer.data_list[cid].double_sided = true;
    viewer.data_list[cid].set_face_based(true);
    viewer.data_list[cid].show_lines = true;
    viewer.data_list[cid].show_faces = true;
    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string matcap_filepath = matcap_file;
    igl::png::readPNG(matcap_filepath, R, G, B, A);
    viewer.data_list[cid].set_face_based(true);
    viewer.data_list[cid].use_matcap = true;
    viewer.data_list[cid].set_texture(R, G, B, A);
    viewer.data_list[cid].show_texture = true;
}

void InteractiveCDHook::set_viewer_clusters(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters, int cid, int fid)
{

    viewer.data_list[fid].clear();
    viewer.data_list[cid ].clear();
    viewer.data_list[cid ].set_mesh(V, F); //should only send V_ext
    viewer.data_list[fid].set_mesh(V_high_res, F_high_res); //should only send V_ext
    viewer.data_list[fid].is_visible = false; //dont wanna see or hear about you when drawing clusters young man!
    viewer.data_list[cid].is_visible = true; //dont wanna see or hear about you when using matcap young man... though tbh why not??

    viewer.data_list[cid].show_face_labels = false;
    viewer.data_list[cid].show_vertex_labels = false;

    viewer.data_list[fid].show_face_labels = false;
    viewer.data_list[fid].show_vertex_labels = false;

    viewer.data_list[cid ].point_size = 10;
          
    viewer.data_list[fid].show_lines = true;
    viewer.data_list[cid ].show_faces = true;
    viewer.data_list[cid ].invert_normals = true;
    viewer.data_list[cid ].double_sided = true;
    viewer.data_list[cid ].set_face_based(true);

    viewer.data_list[cid ].use_matcap = false;
    Eigen::MatrixXd colormap = get_rainbow_colormap();
    viewer.data_list[cid ].set_data(clusters.cast<double>());
    viewer.data_list[cid ].set_colormap(colormap);
}

void InteractiveCDHook::draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{
    refresh = ImGui::Button("Refresh"); 
    if (refresh)
    {
        std::cout << "updating sim state..." << std::endl;
    }
    
   // double sim_fps = 1.0 / (timings.mean());
   
    char string_sim_fps[1024];
    sprintf(string_sim_fps, "Timestep cost (s) : % .3g", timings.mean());   
    ImGui::Text(string_sim_fps);

    ImGui::Text("Vertices : %5d", V.rows());
    ImGui::Text("Tets : %5d", T.rows());

    ImGui::Text("Display Vertices : %5d", V_high_res.rows());
    ImGui::Text("Display Faces : %5d", F_high_res.rows());

    ImGui::Text("Steps : %i", step);
    if (ImGui::CollapsingHeader("Visualization"))
    {
        ImGui::Checkbox("vis CD", &v_state.vis_cd);

        bool use_gpu_proj = as.proj_gpu > 0;
        bool changed_proj = ImGui::Checkbox("Project on GPU: ", &use_gpu_proj);
        if (changed_proj)
        {
            as.proj_gpu = use_gpu_proj ? true : false;
                                    //cant run gpu projection at the same time as any fancy libigl coloring
            if (as.proj_gpu)
            {
                set_viewer_defo_textures(*viewer, V0, F, cd_sim.B, rig->W, v_state.coarse_vis_id, v_state.fine_vis_id);
                
            }
        }
        if (as.proj_gpu == 0)
        {
            as.proj_gpu = 0;

            ImGui::RadioButton("Matcap", (int*)&new_v_state.vis_mode, (int)VIS_MODE::MATCAP); ImGui::SameLine();
            ImGui::RadioButton("Clusters", (int*)&new_v_state.vis_mode, (int)VIS_MODE::CLUSTERS); ImGui::SameLine();
            ImGui::RadioButton("Textures", (int*)&new_v_state.vis_mode, (int)VIS_MODE::TEXTURES);
            if (v_state.vis_mode != new_v_state.vis_mode || changed_proj)
            {
                v_state.vis_mode = new_v_state.vis_mode;
                init_viewer(*viewer);

            }
            if (v_state.vis_mode == TEXTURES)
            {
                ImGui::Checkbox("Draw Cage", &v_state.show_cage);
                viewer->data_list[v_state.coarse_vis_id].is_visible = v_state.show_cage;

            }
        }
    }

    if (ImGui::CollapsingHeader("Rigs"))
    {
       
        if (ImGui::BeginListBox("Saved Rigs"))
        {
            for (int n = 0; n < rig_names.size(); n++)
            {
                const bool is_selected = (as.rig_id == n);
                if (ImGui::Selectable(rig_names[n].c_str(), is_selected))
                    new_as.rig_id = n;

                // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndListBox();
        }
      /*  if (ImGui::BeginListBox("Create Rigs"))
        {

            for (int n = 0; n < create_rig_names.size(); n++)
            {
                const bool is_selected = (as.rig_id == n);
                if (ImGui::Selectable(create_rig_names[n].c_str(), is_selected))
                {
                    if (create_rig_names[n] == "two_handle_rig")
                    {
                        create_two_handle_rig(create_rig_paths[n], V0);
                    }
                }

                // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                if (is_selected)
                    ImGui::SetItemDefaultFocus();

                ImGui::EndListBox();
            }
        }*/
       
    }


    ImGui::Text("Animation Mode");
    ImGui::RadioButton("Play Around", (int*)&as.animation_mode, (int)ANIMATION_MODE::INTERACTIVE_ANIMATION); ImGui::SameLine();
    if (ImGui::RadioButton("Mode Anim", (int*)&as.animation_mode, (int)ANIMATION_MODE::EIGENMODES_ANIMATION))
        step = 0;
    if (as.animation_mode == ANIMATION_MODE::EIGENMODES_ANIMATION)
    {

        std::string current_mode_str = "Current Mode : " + std::to_string(mas.mode);
        ImGui::Text(current_mode_str.c_str());
    }
    bool changed_constraint = ImGui::RadioButton("Pinning", (int*) &as.constraint_type, (int)CONSTRAINT_TYPE::PINNING); ImGui::SameLine();
    changed_constraint = changed_constraint || ImGui::RadioButton("CD", (int*) &as.constraint_type, (int)CONSTRAINT_TYPE::COMPLEMENTARY_DYNAMICS);
    
    if (changed_constraint)
    {
        init_simulation();
        init_viewer(*viewer);
    }
  

    ImGui::SliderFloat("Young's Modulus", &new_as.ym, 0.1, 1000, "% .3f", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Poisson Ratio", &new_as.pr, 0.0, 0.5);
    ImGui::Checkbox("Do Reduction", &new_as.do_reduction); ImGui::SameLine();
    ImGui::Checkbox("Do Clustering", &new_as.do_clustering);
    // ImGui::DragInt("# Modes", &new_num_modes, 1, 500);
     //ImGui::DragInt("# Clusters", &new_num_clusters, 1, 500);
    ImGui::InputInt("# Modes", &new_as.r, 25, 100);
    new_as.r = std::fmax(new_as.r, 1);
    ImGui::InputInt("# Clusters", &new_as.l, 25, 100);
    ImGui::InputInt("# features", &new_as.feat, 1, 10);
    new_as.l = std::fmax(new_as.l, 1);

    if (ImGui::CollapsingHeader("Record"))
    {
        if (ImGui::Button("Save"))
        {
            save_results();
        }
        if (ImGui::Checkbox("metrics", &as.record_metrics))
        {
            time_energy.resize(0, 2);
        }
    }
    as.rig_controller->draw_gui(menu);

}

bool InteractiveCDHook::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {

    return as.rig_controller->mouse_down(viewer, button, modifier);
  //  return rig->mouse_down(viewer, button, modifier);

 //   return false;

}

bool InteractiveCDHook::mouse_move(igl::opengl::glfw::Viewer& viewer, int x, int y) {

    return rig->mouse_move(viewer, x, y);

}

bool InteractiveCDHook::mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    //should call a rig class wrapper
    return rig->mouse_up(viewer, button, modifier);
}

bool InteractiveCDHook::key_callback(igl::opengl::glfw::Viewer& viewer, unsigned int button, int modifier)
{
    if (as.rig_controller->key_callback(viewer, button, modifier))
        return true;
   // if (as.constraint_controller->key_callback(viewer, button, modifier))
   //     return true;
     if (!rig->key_callback(viewer, button, modifier))
    {
        //rewrite this because the default mouse controls dont work with libigl new data_id object
        switch (button)
        {
        case 'F':
        case 'f':
            //why does this cause such a problem, whereas before it didnt?
           // data().set_face_based(!data().face_based);
            viewer.data_list[v_state.coarse_vis_id].set_face_based(!viewer.data_list[v_state.coarse_vis_id].face_based);
            return true;
        case 'T':
        case 't':
            viewer.data_list[v_state.coarse_vis_id].show_faces = !viewer.data_list[v_state.coarse_vis_id].show_faces;
            return true;

        case 'D':
        case 'd':
        {
            viewer.data_list[v_state.coarse_vis_id].double_sided = !viewer.data_list[v_state.coarse_vis_id].double_sided;
            return true;
        }
        case 'I':
        case 'i':
        {
            viewer.data_list[v_state.coarse_vis_id].dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
            viewer.data_list[v_state.coarse_vis_id].invert_normals = !viewer.data_list[v_state.coarse_vis_id].invert_normals;
            return true;
        }
        case 'L':
        case 'l':
        {
            viewer.data_list[v_state.coarse_vis_id].show_lines = !viewer.data_list[v_state.coarse_vis_id].show_lines;
            return true;
        }
        case 's': //enter key is pressed
        {
            refresh = true;
            std::cout << "updating sim state..." << std::endl;
            return true;
        }
        default:
            return false;
        }
    }
    else
        return true;
}


void InteractiveCDHook::change_rig_type()
{
    as.rig_id = new_as.rig_id;

    igl::opengl::glfw::imgui::ImGuizmoWidget* imguizmo = rig->gizmoPlugin;

    init_rig(rig_paths[as.rig_id], as.mesh_file_path);

    as.rig_controller->init_guizmo_viewer(viewer, guizmo);

    init_simulation();


    cd_sim.update_compelementary_constraint(rig->J, as.cd_mode_dir, as.cd_clusters_dir);
    
    init_rig_controller(rig);
    sim.update_equality_constraint(rig->S);

    igl::slice(cd_sim.B, iext, 1, cd_B_ext);
    igl::slice(sim.B, iext, 1, pinned_B_ext);
    igl::slice(cd_sim.J, iext, 1, cd_J_ext);		//slice J the same as B
    init_viewer(*viewer);
}

void InteractiveCDHook::change_animation_mode()
{
    as.animation_mode = new_as.animation_mode;
    init_simulation();

}

/*
Poll all simulation changes if refresh is clicked. This is to avoid recalculating modes as you're moving the slider.
*/
void InteractiveCDHook::poll_sim_changes()
{
    std::cout << "Changing Simulation State..." << std::endl;
    if (as.do_clustering != new_as.do_clustering)
    {
        as.do_clustering = new_as.do_clustering;
        cd_sim.switch_clustering(as.do_clustering);
        sim.switch_clustering(as.do_clustering);
        init_viewer(*viewer);
    }
    if (as.do_reduction != new_as.do_reduction)
    {
        as.do_reduction = new_as.do_reduction;
        cd_sim.switch_reduction(as.do_reduction);
        sim.switch_reduction(as.do_reduction); //never do reduction for pininng method this
        init_viewer(*viewer);
    }
    if (new_as.ym != as.ym || new_as.pr != as.pr)
    {
        as.pr = new_as.pr;
        as.ym = new_as.ym;
        std::cout << "Changed material parameters... recomputing matrices and factorization.." << std::endl;
        cd_sim.update_material_properties(as.ym, as.pr);
        sim.update_material_properties(as.ym, as.pr);
    }

    if (new_as.l != as.l || new_as.feat != as.feat)
    {
        as.l = new_as.l;
        as.feat = new_as.feat;
        sim.num_modal_features = as.feat;

        cd_sim.update_clusters(new_as.l); //TODO add as.feat to cd_sim class
        sim.update_clusters(new_as.l);
        init_viewer(*viewer);
    }
    if (new_as.r != as.r)
    {
        as.r = new_as.r;
        std::cout << "Changed #Modes to " + std::to_string(as.r) + " checking cache..." << std::endl;
        z_curr = Eigen::VectorXd::Zero(as.r);
        z_prev = Eigen::VectorXd::Zero(as.r);
        cd_sim.update_modes(new_as.r);
        sim.update_modes(new_as.r);
        igl::slice(cd_sim.B, iext, 1, cd_B_ext);
        igl::slice(sim.B, iext, 1, pinned_B_ext);
        igl::slice(cd_sim.J, iext, 1, cd_J_ext);		//slice J the same as B
        init_viewer(*viewer);
      
    }
    if (as.rig_id != new_as.rig_id)
    {
        change_rig_type();
    }
    if (new_as.animation_mode != as.animation_mode)
    {
        change_animation_mode();
    }


    std::cout << "Done!" << std::endl;
}