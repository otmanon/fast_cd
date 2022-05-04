#include "FastCDViewer.h"
#include "launch_viewer_custom_shader.h"
#include "rainbow_cmap.h"

#include <igl/LinSpaced.h>
#include <igl/opengl/gl.h>
#include <igl/png/readPNG.h>
FastCDViewer::FastCDViewer()
{
	viewer = new igl::opengl::glfw::Viewer();

	while (viewer->data_list.size() < 3) viewer->append_mesh();

    viewer->core().background_color.setOnes();
	viewer->core().animation_max_fps = 240;
	viewer->core().is_animating = true;
    fid = 1;
    cid = 0;
    viewer->data_list[2].is_visible = true;
}


FastCDViewer::FastCDViewer(igl::opengl::glfw::Viewer* viewer, igl::opengl::glfw::imgui::ImGuiWidget* guizmo)
{
    this->viewer = viewer;
    while (viewer->data_list.size() < 3) viewer->append_mesh();
    viewer->core().animation_max_fps = 240;
    viewer->core().is_animating = true;
    fid = 1;
    cid = 0;
    this->guizmo = guizmo;
}


FastCDViewer::FastCDViewer(igl::opengl::glfw::Viewer* viewer)
{
    this->viewer = viewer;
	while (viewer->data_list.size() < 3) viewer->append_mesh();
    viewer->core().animation_max_fps = 240;
    viewer->core().is_animating = true;
    fid = 1;
    cid = 0;
}
void FastCDViewer::launch()
{
	launch_viewer_custom_shader(*viewer, true, false, "fast CD app", 1920, 1080);
}


void FastCDViewer::configure_clusters(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& clusters)
{
    viewer->data_list[fid].clear();
    viewer->data_list[cid].clear();
    viewer->data_list[cid].set_mesh(V, F); //should only send V_ext
  //  viewer->data_list[fid].set_mesh(V_high_res, F_high_res); //should only send V_ext
    viewer->data_list[fid].is_visible = false; //dont wanna see or hear about you when drawing clusters young man!
    viewer->data_list[cid].is_visible = true; //dont wanna see or hear about you when using matcap young man... though tbh why not??
    //viewer->data_list[2].is_visible = false;
    viewer->data_list[cid].show_face_labels = false;
    viewer->data_list[cid].show_vertex_labels = false;

    viewer->data_list[fid].show_face_labels = false;
    viewer->data_list[fid].show_vertex_labels = false;

    viewer->data_list[cid].point_size = 10;

    viewer->data_list[fid].show_lines = true;
    viewer->data_list[cid].show_faces = true;
    viewer->data_list[cid].invert_normals = true;
    viewer->data_list[cid].double_sided = true;
    viewer->data_list[cid].set_face_based(true);

    viewer->data_list[cid].use_matcap = false;
    Eigen::MatrixXd colormap = get_rainbow_colormap();
    viewer->data_list[cid].set_data(clusters.cast<double>());
    viewer->data_list[cid].set_colormap(colormap);
}

void FastCDViewer::configure_deformation_texture( Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& B, Eigen::MatrixXd& W)
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
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> tex = Eigen::Matrix< float, Eigen::Dynamic, 3, Eigen::RowMajor>::Zero(s * s, 3);

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
   // igl::opengl::glfw::Viewer viewer;

    viewer->data_list[cid].clear();
    viewer->data_list[fid].clear();
    viewer->data_list[cid].set_mesh(V, F);
    viewer->data_list[cid].set_face_based(false);
    viewer->data_list[cid].use_matcap = false;
    viewer->data_list[cid].show_lines = true;
    viewer->data_list[cid].show_faces = true;
    viewer->data_list[cid].invert_normals = false;
    viewer->data_list[cid].double_sided = true;
    viewer->data_list[cid].show_lines = false;
    viewer->data_list[cid].show_texture = false;
    viewer->data_list[fid].show_texture = false;
    viewer->data_list[fid].is_visible = false;
    viewer->data_list[cid].is_visible = true;


    viewer->data_list[cid].show_face_labels = false;
    viewer->data_list[cid].show_vertex_labels = false;

    viewer->data_list[fid].show_face_labels = false;
    viewer->data_list[fid].show_vertex_labels = false;
    ///////////////////////////////////////////////////////////////////
   // Send texture and vertex attributes to GPU... should make this its own function but maybe it's okay
   ///////////////////////////////////////////////////////////////////
    {
        GLuint prog_id = viewer->data_list[cid].meshgl.shader_mesh;
        glUseProgram(prog_id);
        GLuint VAO = viewer->data_list[cid].meshgl.vao_mesh;
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
        glBindTexture(GL_TEXTURE_2D, viewer->data_list[cid].meshgl.vbo_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, s, s, 0, GL_RGB, GL_FLOAT, tex.data());
    }
    //This is very important... otherwise when the new mesh is set, it triggers a dirty call , and igl overwrites the texture with junk...
    viewer->data_list[cid].dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;



}

void FastCDViewer::configure_color_texture(std::string texture_filepath, Eigen::MatrixXd& V_coarse, Eigen::MatrixXi& F_coarse, Eigen::MatrixXd& V_fine, Eigen::MatrixXi& F_fine,
    Eigen::MatrixXd& UV_fine, Eigen::MatrixXi& FUV_fine)
{

    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh
    viewer->data_list[fid].clear();
    viewer->data_list[cid].clear();
   // viewer->data_list[2].is_visible = true;
    viewer->data_list[cid].show_lines = false;
    viewer->data_list[cid].show_faces = false;
    viewer->data_list[cid].is_visible = false;//dont wanna see or hear about you when using textures
    viewer->data_list[cid].show_face_labels = false;
    viewer->data_list[cid].show_vertex_labels = false;
    viewer->data_list[cid].set_mesh(V_coarse, F_coarse);



    viewer->data_list[fid].is_visible = true;
    viewer->data_list[fid].invert_normals = false;
    viewer->data_list[fid].double_sided = false;
    viewer->data_list[fid].set_face_based(true);
    viewer->data_list[fid].show_vertex_labels = false;
    viewer->data_list[fid].show_face_labels = false;
    viewer->data_list[fid].set_mesh(V_fine, F_fine);

    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string texture_path = texture_filepath;
    bool read = igl::png::readPNG(texture_path, R, G, B, A);

    if (UV_fine.rows() > 0)  //if our texture matches the display mesh
    {
        //  viewer->data().set_normals(FN);
        viewer->data_list[fid].set_uv(UV_fine, FUV_fine);
        viewer->data_list[fid].show_lines = false;
        viewer->data_list[fid].show_texture = true;
        viewer->data_list[fid].show_faces = true;
        viewer->data_list[fid].set_texture(R, G, B);
    }
    else
    {
        viewer->data_list[fid].invert_normals = true;
    }

}


void FastCDViewer::configure_matcap(std::string matcap_file, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{

    viewer->data_list[fid].clear();
    viewer->data_list[cid].clear();

    viewer->data_list[cid].set_mesh(V, F);
    //viewer->data_list[fid].set_mesh(V_high_res, F_high_res); //should only send V_ext

    viewer->data_list[fid].is_visible = false; //dont wanna see or hear about you when using matcap young man... though tbh why not??
    viewer->data_list[cid].is_visible = true; //dont wanna see or hear about you when using matcap young man... though tbh why not??


    viewer->data_list[cid].point_size = 10;


    viewer->data_list[cid].show_face_labels = false;
    viewer->data_list[cid].show_vertex_labels = false;

    viewer->data_list[fid].show_face_labels = false;
    viewer->data_list[fid].show_vertex_labels = false;

    viewer->data_list[cid].invert_normals = false;
    viewer->data_list[cid].double_sided = true;
    viewer->data_list[cid].set_face_based(true);
    viewer->data_list[cid].show_lines = true;
    viewer->data_list[cid].show_faces = true;
    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string matcap_filepath = matcap_file;
    igl::png::readPNG(matcap_filepath, R, G, B, A);
    viewer->data_list[cid].set_face_based(true);
    viewer->data_list[cid].use_matcap = true;
    viewer->data_list[cid].set_texture(R, G, B, A);
    viewer->data_list[cid].show_texture = true;
}


void FastCDViewer::render_full(Eigen::MatrixXd& V, int id)
{
    GLuint prog_id = viewer->data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, false); //dopn't use gpu to compute displacements!!
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);

    viewer->data_list[id].set_vertices(V);

}

void FastCDViewer::render_reduced_cpu_proj(Eigen::VectorXd& z, Eigen::VectorXd& p, Eigen::MatrixXd& B, Eigen::SparseMatrix<double>& J, int id)
{
    GLuint prog_id = viewer->data_list[cid].meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, false); //dopn't use gpu to compute displacements!!
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);
    Eigen::VectorXd u = B * z.cast<double>() + J * p.cast<double>();  //the titular cpu proj
    Eigen::MatrixXd V = Eigen::Map<Eigen::MatrixXd>(u.data(), u.rows() / 3, 3);

    viewer->data_list[id].set_vertices(V);
}


void FastCDViewer::render_reduced_gpu_proj(Eigen::VectorXd& z, Eigen::VectorXd& p, int n, int cid)
{
    //TODO ensure that the viewer is configured for texture
      //TODO: assert that viewer is already configured for reduced step... otherwise need to resend texture
    Eigen::VectorXf zf = z.cast<float>();
    Eigen::VectorXf pf = p.cast<float>();
    const int m = z.rows();
    const int b = int(p.rows() / 12);
    const int s = ceil(sqrt(n * (m + b)));

    GLuint prog_id = viewer->data_list[cid].meshgl.shader_mesh;
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, viewer->data_list[cid].meshgl.vbo_tex);
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
    glUniform1fv(q_loc, m, zf.data());
    GLint p_loc = glGetUniformLocation(prog_id, "p"); //rig parameters
    glUniform1fv(p_loc, p.rows(), pf.data());
    GLint cd_loc = glGetUniformLocation(prog_id, "proj_gpu");
    glUniform1i(cd_loc, true); //if we call this then we must assume it's true
    // Do this now so that we can stop texture from being loaded by viewer
    GLint pin_loc = glGetUniformLocation(prog_id, "pin");
    glUniform1i(pin_loc, false);
    if (viewer->data_list[cid].dirty)
    {
        viewer->data_list[cid].updateGL(
            viewer->data_list[cid],
            viewer->data_list[cid].invert_normals,
            viewer->data_list[cid].meshgl
        );
        viewer->data_list[cid].dirty = igl::opengl::MeshGL::DIRTY_NONE;
    }

    viewer->data_list[cid].dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;

}

void draw_gui(igl::opengl::glfw::imgui::ImGuiMenu& menu);

void FastCDViewer::set_pre_draw_callback(std::function<bool()>& callback)
{
	viewer->callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
	{
		callback();
		return false;
	};
}
