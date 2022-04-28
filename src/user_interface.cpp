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

#include "igl/opengl/gl.h"
#include "igl/LinSpaced.h"
#include "igl/opengl/create_shader_program.h"
#include "igl/read_triangle_mesh.h"
#include "igl/readDMAT.h"
#include "igl/readPLY.h"
//Render UI related
void InteractiveCDHook::render(igl::opengl::glfw::Viewer& viewer)
{
  //  Eigen::VectorXf q0 = Eigen::VectorXf::Zero(as.r);
    Eigen::VectorXf z = 0*z_curr.cast<float>();
    z(0) = 1;
  //  const int keyrate = 10;
   // if (step % keyrate == 0)
  //  {
       // q0 = q1;
       // z_curr(0) += 1;// = z.cast<double>();
     //   q1 = q1.array().pow(100.0).eval();
   //}
    
  //  Eigen::VectorXf qa = q0 + double(count % keyrate) / (keyrate - 1.0) * (q1 - q0);
   // qa /= qa.sum();
   // count++;
    /////////////////////////////////////////////////////////
    // Send uniforms to shader
    /////////////////////////////////////////////////////////
    const int s = ceil(sqrt(V.rows() * as.r));

    GLuint prog_id = viewer.data().meshgl.shader_mesh;
    glUseProgram(prog_id);
    GLint n_loc = glGetUniformLocation(prog_id, "n");
    glUniform1i(n_loc, V.rows());
    GLint m_loc = glGetUniformLocation(prog_id, "m");
    glUniform1i(m_loc, as.r);
    GLint s_loc = glGetUniformLocation(prog_id, "s");
    glUniform1i(s_loc, s);
    GLint q_loc = glGetUniformLocation(prog_id, "q");
    glUniform1fv(q_loc, sim.B.cols(), z.data());
    // Do this now so that we can stop texture from being loaded by viewer
   
        viewer.data().updateGL(viewer.data(), viewer.data().invert_normals, viewer.data().meshgl);
        viewer.data().dirty = igl::opengl::MeshGL::DIRTY_NONE;
        viewer.data().dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;
   // viewer.data().dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;

   // as.rig_controller->render(viewer);
   // rig->render(viewer);
    
   // if (v_state.vis_mode == TEXTURES)
   // {
   //     viewer.data_list[v_state.fine_vis_id].set_vertices(V_high_res);
   // }
   // 
   // viewer.data_list[v_state.coarse_vis_id].set_vertices(V_ext);
   

   //viewer.data_list[v_state.coarse_vis_id].set_vertices(V_ext);
   // viewer.data_list[v_state.coarse_vis_id].set_vertices(V_high_res);

    //maybe can let the fragment shader take care of this
  //  viewer.data_list[v_state.coarse_vis_id].compute_normals();d

   
}

void InteractiveCDHook::init_viewer(igl::opengl::glfw::Viewer& bogusViewer)
{
    using namespace Eigen;
    using namespace std;
  //  Eigen::MatrixXd V;
   // Eigen::MatrixXi F;
   // igl::readPLY("C:/Users/otman/OneDrive/Desktop/matrixdefo/face.ply", V, F);

    ///////////////////////////////////////////////////////////////////
    // Load and prepare data
    ///////////////////////////////////////////////////////////////////
    Eigen::Matrix< float, Eigen::Dynamic, 1> I;
    Eigen::Matrix< float, Eigen::Dynamic, 3, Eigen::RowMajor> tex;
   // Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> U;
   // U = Eigen::Map<Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(sim.B.data(), )
   // igl::readDMAT("C:/Users/otman/OneDrive/Desktop/matrixdefo/data.DMAT", U);
    assert((cd_sim.B.rows() == V.rows() * 3) && "#U should be 3*#V");
    //std::cout<<"**warning** resizing to min(U.cols(),100)"<<std::endl;
    //U.conservativeResize(U.rows(),std::min(100,(int)U.cols()));
    I = igl::LinSpaced< Eigen::Matrix< float, Eigen::Dynamic, 1> >(V.rows(), 0, V.rows() - 1);
    const int n = V.rows();
    const int m = cd_sim.B.cols();
    const int s = ceil(sqrt(n)); // ceil(sqrt(n * m));
   // assert(s * s > n * m);//
    printf("%d %d %d\n", n, m, s);
    tex = Eigen::MatrixXf::Zero(n, 3);
    Eigen::VectorXd x = V.col(0);
    x.array() /= (x.maxCoeff() - x.minCoeff());
    x.array() -= x.minCoeff();
    tex.block(0, 0, V.rows(), 1) = x.cast<float>(); //first column is all x.
   //for (int j = 0; j < m; j++)
   //{
   //    for (int i = 0; i < n; i++)
   //    {
   //        for (int c = 0; c < 3; c++)
   //        {
   //            tex(m*i + j, c) = -float(cd_sim.B(i + c * n, j)); // should give zero displacement for all modes
   //        }
   //    }
   //}

     ///////////////////////////////////////////////////////////////////
    // Compile Shaders
    ///////////////////////////////////////////////////////////////////
    std::string mesh_vertex_shader_string =
        R"(#version 150
      uniform mat4 view;
      uniform mat4 proj;
      uniform mat4 normal_matrix;
      in vec3 position;
      in vec3 normal;
      out vec3 position_eye;
      out vec3 normal_eye;
      out vec4 tex_color;
      in vec4 Ka;
      in vec4 Kd;
      in vec4 Ks;
     // in vec2 texcoord;
      out vec2 texcoordi;
      out vec4 Kai;
      out vec4 Kdi;
      out vec4 Ksi;

     in float id;
     uniform int n;
     uniform int m;
     uniform int s;
     uniform float q[512];
    uniform sampler2D tex;

    void main()
    {
       vec3 displacement = vec3(0,0,0);
    // for(int j = 0;j < m; j++)
    // {
    //   int index = int(id)*m+j;
    //   int si = index % s;
    //   int sj = int((index - si)/s);
    //   displacement = displacement + texelFetch(tex,ivec2(si,sj),0).xyz*q[j];
    // }
      int index= int(id);
      int si = index % s;
      int sj = int((index - si)/s);
      vec3 deformed = position;// + displacement ; 
      tex_color = vec4(texelFetch(tex,ivec2(si,sj),0).xyz,  1.0);
      position_eye =  vec3 (view * vec4 (deformed, 1.0));
      gl_Position = proj * vec4 (position_eye, 1.0);
      normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
      normal_eye = normalize(normal_eye);
      Kai = Ka;
      Kdi = Kd;
      Ksi = Ks;
    })";


    std::string mesh_fragment_shader_string =
        R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform vec4 fixed_color;
  in vec3 position_eye;
  in vec3 normal_eye;
  uniform vec3 light_position_eye;
  vec3 Ls = vec3 (1, 1, 1);
  vec3 Ld = vec3 (1, 1, 1);
  vec3 La = vec3 (1, 1, 1);
  in vec4 Ksi;
  in vec4 Kdi;
  in vec4 Kai;
  in vec4 tex_color;
  in vec2 texcoordi;
  uniform sampler2D tex;
  uniform float specular_exponent;
  uniform float lighting_factor;
  uniform float texture_factor;
  uniform float matcap_factor;
  uniform float double_sided;
  out vec4 outColor;
  void main()
  {

    vec3 xTangent = dFdx( position_eye );
    vec3 yTangent = dFdy( position_eye );
    vec3 faceNormal = normalize( cross( xTangent, yTangent ) );
    if(matcap_factor == 1.0f)
    {
      vec2 uv = normalize(faceNormal).xy * 0.5 + 0.5;
      outColor = texture(tex, uv);
    }else
    {
      vec3 Ia = La * vec3(Kai);    // ambient intensity

      vec3 vector_to_light_eye = light_position_eye - position_eye;
      vec3 direction_to_light_eye = normalize (vector_to_light_eye);
      float dot_prod = dot (direction_to_light_eye, normalize(faceNormal));
      float clamped_dot_prod = abs(max (dot_prod, -double_sided));
      vec3 Id = Ld * vec3(Kdi) * clamped_dot_prod;    // Diffuse intensity

      vec3 reflection_eye = reflect (-direction_to_light_eye, normalize(faceNormal));
      vec3 surface_to_viewer_eye = normalize (-position_eye);
      float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
      dot_prod_specular = float(abs(dot_prod)==dot_prod) * abs(max (dot_prod_specular, -double_sided));
      float specular_factor = pow (dot_prod_specular, specular_exponent);
      vec3 Is = Ls * vec3(Ksi) * specular_factor;    // specular intensity
      vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);
      outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;
      
      if (fixed_color != vec4(0.0)) outColor = fixed_color;
     
     outColor = tex_color;
    }
  }
)";


    ///////////////////////////////////////////////////////////////////
    // Initialize viewer and opengl context
    ///////////////////////////////////////////////////////////////////
    Eigen::MatrixXf texT = tex.transpose();
    igl::opengl::glfw::Viewer v;
    v.launch_init(true, false, "fast CD App", 1920, 1080);
    //destroys any existing shader programs
    v.data().meshgl.free();
    v.data().meshgl.is_initialized = true;
    v.data().meshgl.init_buffers();
  //  v.data().meshgl.init_text_rendering();
    igl::opengl::create_shader_program(
        mesh_vertex_shader_string,
        mesh_fragment_shader_string,
        {},
        v.data().meshgl.shader_mesh);
    ///////////////////////////////////////////////////////////////////
    // Send texture and vertex attributes to GPU
    ///////////////////////////////////////////////////////////////////
    {
        GLuint prog_id = v.data().meshgl.shader_mesh;
        glUseProgram(prog_id);
        GLuint VAO = v.data().meshgl.vao_mesh;
        glBindVertexArray(VAO);
        GLuint IBO;
        glGenBuffers(1, &IBO);
        glBindBuffer(GL_ARRAY_BUFFER, IBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * I.size(), I.data(), GL_STATIC_DRAW);
        GLint iid = glGetAttribLocation(prog_id, "id");
        glVertexAttribPointer(
            iid, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(iid);
        glBindVertexArray(0);
        glActiveTexture(GL_TEXTURE0);
        //glGenTextures(1, &v.opengl.vbo_tex);
        glBindTexture(GL_TEXTURE_2D, v.data().meshgl.vbo_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        // 8650×8650 texture was roughly the max I could still get 60 fps, 8700²
        // already dropped to 1fps
        //
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, s, s, 0, GL_RGB, GL_FLOAT, tex.data());
    }


    Eigen::VectorXd b = cd_sim.B.col(0);
    Eigen::MatrixXd UB = Eigen::Map<Eigen::MatrixXd>(b.data(), V.rows(), 3);
    Eigen::MatrixXd U = 0.01*UB + V0;
    //  v.data().set_mesh(U, F);
    v.data().set_mesh(V0, F);

    Eigen::VectorXf q0 = Eigen::VectorXf::Zero(m, 1);
    q0(0) = 0.1;
    Eigen::VectorXf q1 = Eigen::VectorXf::Zero(m, 1);

    int i = 0;
    v.callback_pre_draw = [&](igl::opengl::glfw::Viewer& v) ->bool
    {
       //static size_t count = 0;
       //const int keyrate = 15;
       //if (count % keyrate == 0)
       //{
       //    q0 = q1;
       //    q1 = Eigen::VectorXf::Random(m, 1).array() * 0.5 + 0.5;
       //    q1 = q1.array().pow(100.0).eval();
       //}
        //Eigen::VectorXf qa = q0 + double(count % keyrate) / (keyrate - 1.0) * (q1 - q0);
        //qa /= qa.sum();
        //count++;
        /////////////////////////////////////////////////////////
        // Send uniforms to shader
        /////////////////////////////////////////////////////////
        GLuint prog_id = v.data().meshgl.shader_mesh;
        glUseProgram(prog_id);
        GLint n_loc = glGetUniformLocation(prog_id, "n");
        glUniform1i(n_loc, n);
        GLint m_loc = glGetUniformLocation(prog_id, "m");
        glUniform1i(m_loc, m);
        GLint s_loc = glGetUniformLocation(prog_id, "s");
        glUniform1i(s_loc, s);
        GLint q_loc = glGetUniformLocation(prog_id, "q");
        glUniform1fv(q_loc, cd_sim.B.cols(), q0.data());
       // Eigen::VectorXd b = cd_sim.B.col(0);
       // Eigen::MatrixXd UB = Eigen::Map<Eigen::MatrixXd>(b.data(), V.rows(), 3);
       // Eigen::MatrixXd U = UB + V0;
      //  v.data().set_mesh(U, F);
        // Do this now so that we can stop texture from being loaded by viewer
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, s, s, 0, GL_RGB, GL_FLOAT, tex.data());
        if (v.data().dirty)
        {
            v.data().updateGL(v.data(), v.data().invert_normals, v.data().meshgl);
            v.data().dirty = igl::opengl::MeshGL::DIRTY_NONE;
           
        }
        v.data().dirty &= ~igl::opengl::MeshGL::DIRTY_TEXTURE;
        i += 1;
        return false;
    };
    //v.data().show_texture = true;
    v.core().animation_max_fps = 240;
    v.core().is_animating = true;
    v.launch_rendering(true);
    v.launch_shut();
    
}


void InteractiveCDHook::set_viewer_textures()
{
    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh
    viewer->data_list[v_state.fine_vis_id].clear();
    viewer->data_list[v_state.coarse_vis_id].show_lines = true;
    viewer->data_list[v_state.coarse_vis_id].show_faces = false;
    viewer->data_list[v_state.coarse_vis_id].set_mesh(V_ext, F_ext);

    viewer->data_list[v_state.fine_vis_id].clear();
        

    viewer->data_list[v_state.fine_vis_id].invert_normals = false;
    viewer->data_list[v_state.fine_vis_id].double_sided = false;
    viewer->data_list[v_state.fine_vis_id].set_face_based(true);

    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string texture_path = as.texture_png_path;
    bool read = igl::png::readPNG(texture_path, R, G, B, A);

    viewer->data_list[v_state.fine_vis_id].set_mesh(V_high_res, F_high_res);
    if (UV_high_res.rows() > 0)  //if our texture matches the display mesh
    {
        //  viewer.data().set_normals(FN);
        viewer->data_list[v_state.fine_vis_id].set_uv(UV_high_res, FUV_high_res);
        viewer->data_list[v_state.fine_vis_id].show_lines = false;
        viewer->data_list[v_state.fine_vis_id].show_texture = true;
        viewer->data_list[v_state.fine_vis_id].show_faces = true;
        viewer->data_list[v_state.fine_vis_id].set_texture(R, G, B);
    }else
    {
        viewer->data_list[v_state.fine_vis_id].invert_normals = true;
    }
    
}


void InteractiveCDHook::set_viewer_matcap()
{
    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh

    viewer->data_list[v_state.fine_vis_id].clear();
    viewer->data_list[v_state.coarse_vis_id].clear();
 
    viewer->data_list[v_state.coarse_vis_id].set_mesh(V_ext, F_ext);


    viewer->data_list[v_state.coarse_vis_id].point_size = 10;


    viewer->data_list[v_state.coarse_vis_id].invert_normals = true;
    viewer->data_list[v_state.coarse_vis_id].double_sided = false;
    viewer->data_list[v_state.coarse_vis_id].set_face_based(true);
    viewer->data_list[v_state.coarse_vis_id].show_lines = true;
    viewer->data_list[v_state.coarse_vis_id].show_faces = true;
    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string matcap_filepath = "../data/" + as.matcap_file;
    igl::png::readPNG(matcap_filepath, R, G, B, A);
    viewer->data_list[v_state.coarse_vis_id].set_face_based(true);
    viewer->data_list[v_state.coarse_vis_id].use_matcap = true;
    viewer->data_list[v_state.coarse_vis_id].set_texture(R, G, B, A);
    viewer->data_list[v_state.coarse_vis_id].show_texture = true;
}

void InteractiveCDHook::set_viewer_clusters()
{

    viewer->data_list[v_state.fine_vis_id].clear();
    viewer->data_list[v_state.coarse_vis_id].clear();
    viewer->data_list[v_state.coarse_vis_id].set_mesh(V_ext, F_ext);
   

    viewer->data_list[v_state.coarse_vis_id].point_size = 10;

    viewer->data_list[v_state.coarse_vis_id].show_lines = true;
    viewer->data_list[v_state.coarse_vis_id].show_faces = true;
    viewer->data_list[v_state.coarse_vis_id].invert_normals = true;
    viewer->data_list[v_state.coarse_vis_id].double_sided = false;
    viewer->data_list[v_state.coarse_vis_id].set_face_based(true);

    viewer->data_list[v_state.coarse_vis_id].use_matcap = false;
    Eigen::MatrixXd colormap = get_rainbow_colormap();
    Eigen::VectorXi labels_faces;
    igl::slice(sim.labels, FiT, labels_faces);
    viewer->data_list[v_state.coarse_vis_id].set_data(labels_faces.cast<double>());
    viewer->data_list[v_state.coarse_vis_id].set_colormap(colormap);
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
        ImGui::RadioButton("Matcap", (int*)&new_v_state.vis_mode, (int)VIS_MODE::MATCAP); ImGui::SameLine();
        ImGui::RadioButton("Clusters", (int*)&new_v_state.vis_mode, (int)VIS_MODE::CLUSTERS); ImGui::SameLine();
        ImGui::RadioButton("Textures", (int*)&new_v_state.vis_mode, (int)VIS_MODE::TEXTURES);
        if (v_state.vis_mode != new_v_state.vis_mode)
        {
            v_state.vis_mode = new_v_state.vis_mode;
            if (v_state.vis_mode == VIS_MODE::MATCAP)
            {
                set_viewer_matcap();
            }
            else if (v_state.vis_mode == VIS_MODE::CLUSTERS)
            {
                set_viewer_clusters();
            }
            else if (v_state.vis_mode == VIS_MODE::TEXTURES)
            {
                set_viewer_textures();
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
    ImGui::RadioButton("Mode Anim", (int*)&as.animation_mode, (int)ANIMATION_MODE::EIGENMODES_ANIMATION);
    if (new_as.animation_mode == ANIMATION_MODE::EIGENMODES_ANIMATION)
    {
        std::string current_mode_str = "Current Mode : " + std::to_string(mas.mode);
        ImGui::Text(current_mode_str.c_str());
    }
    bool changed_constraint = ImGui::RadioButton("Pinning", (int*) &as.constraint_type, (int)CONSTRAINT_TYPE::PINNING); ImGui::SameLine();
    changed_constraint = changed_constraint || ImGui::RadioButton("CD", (int*) &as.constraint_type, (int)CONSTRAINT_TYPE::COMPLEMENTARY_DYNAMICS);
    
    if (changed_constraint)
    {
        init_simulation();
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
    rig->draw_gui(menu);
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

    render(*viewer);
    
    cd_sim.update_compelementary_constraint(rig->J, as.cd_mode_dir, as.cd_clusters_dir);
    sim.update_equality_constraint(rig->S);

    igl::slice(cd_sim.B, iext, 1, cd_B_ext);
    igl::slice(sim.B, iext, 1, pinned_B_ext);
    igl::slice(cd_sim.J, iext, 1, cd_J_ext);		//slice J the same as B
    if (v_state.vis_mode == VIS_MODE::MATCAP)
    {
        set_viewer_matcap();
    }
    else if (v_state.vis_mode == VIS_MODE::CLUSTERS)
    {
        set_viewer_clusters();
    }
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
        Eigen::VectorXi labels_faces;
        igl::slice(cd_sim.labels, FiT, labels_faces);
        if (v_state.vis_mode == VIS_MODE::CLUSTERS)
        {
            viewer->data_list[v_state.coarse_vis_id].set_data(labels_faces.cast<double>());
        }
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
        Eigen::VectorXi labels_faces;
        igl::slice(sim.labels, FiT, labels_faces);
        if (v_state.vis_mode == VIS_MODE::CLUSTERS)
        {
            viewer->data_list[v_state.coarse_vis_id].set_data(labels_faces.cast<double>());
        }
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