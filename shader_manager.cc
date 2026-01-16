
#include "shader_manager.h"
#include "oculus_sdk4.h" // For ST_OCULUS definition if needed? or just pass integer/bool
#include "utils.h"

using namespace std;

void ShaderManager::clear() {
  fractal.clear();
  effects.clear();
  dof.clear();
  de_shader.clear();
  fxaa.clear();
}

bool ShaderManager::loadFractal(const string &defines,
                                const char *extra_define) {
  string vertex(default_vs);
  string fragment(default_fs);

  readFile(VERTEX_SHADER_FILE, &vertex);
  readFile(FRAGMENT_SHADER_FILE, &fragment);

  // Process synthetic #include statements.
  size_t inc_pos;
  while ((inc_pos = fragment.find("#include ")) != string::npos) {
    size_t name_start = fragment.find("\"", inc_pos) + 1;
    size_t name_end = fragment.find("\"", name_start);
    string name(fragment, name_start, name_end - name_start);
    size_t line_end = fragment.find("\n", inc_pos);
    string inc;
    readFile(name, &inc);
    fragment.replace(inc_pos, line_end - inc_pos, inc);
  }

  if (extra_define == NULL) {
    glsl_source.assign(defines + fragment);
    return fractal.compile(defines, vertex, fragment) != 0;
  } else {
    string extra_defines = extra_define + defines;
    return fractal.compile(extra_defines, vertex, fragment) != 0;
  }
}

bool ShaderManager::loadDEShader(const string &defines,
                                 const char *extra_define) {
  string vertex(default_vs);
  string fragment(default_fs);

  readFile(VERTEX_SHADER_FILE, &vertex);
  readFile(FRAGMENT_SHADER_FILE, &fragment);

  // Process synthetic #include statements.
  size_t inc_pos;
  while ((inc_pos = fragment.find("#include ")) != string::npos) {
    size_t name_start = fragment.find("\"", inc_pos) + 1;
    size_t name_end = fragment.find("\"", name_start);
    string name(fragment, name_start, name_end - name_start);
    size_t line_end = fragment.find("\n", inc_pos);
    string inc;
    readFile(name, &inc);
    fragment.replace(inc_pos, line_end - inc_pos, inc);
  }

  // We don't update glsl_source for DE shader usually? Or do we?
  // Original only updated glsl_source when extra_define was NULL (main
  // fractal).

  if (extra_define == NULL) {
    return de_shader.compile(defines, vertex, fragment) != 0;
  } else {
    string extra_defines = extra_define + defines;
    return de_shader.compile(extra_defines, vertex, fragment) != 0;
  }
}

bool ShaderManager::loadHelpers(const string &defines, int stereoMode) {
  string vertex(effects_default_vs);
  string fragment(effects_default_fs);

  readFile(EFFECTS_VERTEX_SHADER_FILE, &vertex);
  readFile(EFFECTS_FRAGMENT_SHADER_FILE, &fragment);

  glsl_source.append(fragment);

  bool ok = (effects.compile(defines, vertex, fragment) != 0);

#if 0
  // ST_OCULUS is an enum in boxplorer2.cc, hard to use here without full loop dependency or moving enum.
  // Assuming 6 is ST_OCULUS for now, or we should move StereoMode to a header.
  if (stereoMode == 6) { 
    return ok;
  }
#endif

  string fxaa_vertex;
  string fxaa_fragment;

  readFile(FXAA_VERTEX_SHADER_FILE, &fxaa_vertex);
  readFile(FXAA_FRAGMENT_SHADER_FILE, &fxaa_fragment);

  if (!fxaa_vertex.empty() && !fxaa_fragment.empty()) {
    ok &= (fxaa.compile(defines, fxaa_vertex, fxaa_fragment) != 0);
  }

  string dof_vertex;
  string dof_fragment;

  readFile(DOF_VERTEX_SHADER_FILE, &dof_vertex);
  readFile(DOF_FRAGMENT_SHADER_FILE, &dof_fragment);

  if (!dof_vertex.empty() && !dof_fragment.empty()) {
    ok &= (dof.compile(defines, dof_vertex, dof_fragment) != 0);
  }

  return ok;
}
