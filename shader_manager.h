
#pragma once
#include "shader.h"
#include <string>

#define VERTEX_SHADER_FILE "vertex.glsl"
#define FRAGMENT_SHADER_FILE "fragment.glsl"
#define EFFECTS_VERTEX_SHADER_FILE "effects_vertex.glsl"
#define EFFECTS_FRAGMENT_SHADER_FILE "effects_fragment.glsl"
#define DOF_VERTEX_SHADER_FILE "dof_vertex.glsl"
#define DOF_FRAGMENT_SHADER_FILE "dof_fragment.glsl"
#define FXAA_VERTEX_SHADER_FILE "fxaa_vertex.glsl"
#define FXAA_FRAGMENT_SHADER_FILE "fxaa_fragment.glsl"

class ShaderManager {
public:
  Shader fractal;
  Shader effects;
  Shader dof;
  Shader fxaa;
  Shader de_shader;

  std::string glsl_source;

  void clear();
  bool loadFractal(const std::string &defines, const char *extra_define = NULL);
  bool loadDEShader(const std::string &defines,
                    const char *extra_define = NULL);
  bool loadHelpers(const std::string &defines,
                   int stereoMode); // Need stereoMode for Oculus check? maybe
                                    // check internally or pass bool

private:
  const char *default_vs = "void main() { gl_Position = gl_Vertex; }";
  const char *default_fs =
      "void main() { gl_FragColor = vec4(1.0, 0.0, 1.0, 1.0); }";
  const char *effects_default_vs = "void main() { gl_Position = gl_Vertex; "
                                   "gl_TexCoord[0] = gl_MultiTexCoord0; }";
  const char *effects_default_fs =
      "uniform sampler2D tex; void main() { gl_FragColor = texture2D(tex, "
      "gl_TexCoord[0].st); }";
};
