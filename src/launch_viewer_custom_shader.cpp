#include "launch_viewer_custom_shader.h"
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/create_shader_program.h>


void launch_viewer_custom_shader(igl::opengl::glfw::Viewer& v,
    bool resizeable, bool fullscreen, std::string name, int width, int height)
{

    //v.launch_init(true, false);

    std::string mesh_vertex_shader_string =
        R"(#version 150
      uniform mat4 view;
      uniform mat4 proj;
      uniform mat4 normal_matrix;
      in vec3 position;
      in vec3 normal;
      out vec3 position_eye;
      out vec3 normal_eye;
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
        int index = int(id)*m + 0;
        int si = index%s;
        int sj = int((index-si)/s);
         displacement = displacement + texelFetch(tex,ivec2(si,sj),0).xyz*q[0]; //
       // for(int j = 0;j < m; j++)
       // {
       //   int index = int(id)*m+j;
       //   int si = index % s;
       //   int sj = int((index - si)/s);
       //   displacement = displacement + texelFetch(tex,ivec2(si,sj),0).xyz*q[j];
       // }
      vec3 deformed = position + displacement ; 

      position_eye =  vec3 (view * vec4 (deformed, 1.0));
      gl_Position = proj * vec4 (position_eye, 1.0);
      normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
      normal_eye = normalize(normal_eye);
      Kai = Ka;
      Kdi = Kd;
      Ksi = Ks;
    })";

    /*
    std::string mesh_vertex_shader_string =
        R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform mat4 normal_matrix;
  in vec3 position;
  in vec3 normal;
  out vec3 position_eye;
  out vec3 normal_eye;
  in vec4 Ka;
  in vec4 Kd;
  in vec4 Ks;
  in vec2 texcoord;
  out vec2 texcoordi;
  out vec4 Kai;
  out vec4 Kdi;
  out vec4 Ksi;

  void main()
  {
    position_eye = vec3 (view * vec4 (position, 1.0));
    normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
    normal_eye = normalize(normal_eye);
    gl_Position = proj * vec4 (position_eye, 1.0); //proj * view * vec4(position, 1.0);"
    Kai = Ka;
    Kdi = Kd;
    Ksi = Ks;
    texcoordi = texcoord;
  }
)";*/

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
     
  //  outColor = vec4(0,1, 1, 1);
    }
  }
)";

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
 
    v.launch_init(true, false, "fast CD App", 1920, 1080);
    //destroys any existing shader programs
    v.data().meshgl.free();
    v.data().meshgl.is_initialized = true;
    v.data().meshgl.init_buffers();
    v.data().meshgl.init_text_rendering();
    igl::opengl::create_shader_program(
        mesh_vertex_shader_string,
        mesh_fragment_shader_string,
        {},
        v.data().meshgl.shader_mesh);
 //  igl::opengl::create_shader_program(
 //      overlay_vertex_shader_string,
 //      overlay_fragment_shader_string,
 //      {},
 //      v.data().meshgl.shader_overlay_lines);
 //  igl::opengl::create_shader_program(
 //      overlay_vertex_shader_string,
 //      overlay_point_fragment_shader_string,
 //      {},
 //      v.data().meshgl.shader_overlay_points);
 //  igl::opengl::create_shader_program(
 //      text_geom_shader,
 //      text_vert_shader,
 //      text_frag_shader,
 //      {},
 //      v.data().meshgl.shader_text);

    v.core().animation_max_fps = 240;
    v.launch_rendering(true);
    v.launch_shut();

}