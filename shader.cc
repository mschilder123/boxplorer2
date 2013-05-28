#include "shader.h"

#define GL_DECLARE_ONLY
#include "shader_procs.h"

using namespace std;

void Shader::clear() {
  if (glIsProgram(program_)) {
    GLuint shaders[2];
    GLsizei count = 2;   
    glGetAttachedShaders(program_, count, &count, shaders);
    for (GLsizei i = 0; i < count; ++i) {
      glDetachShader(program_, shaders[i]);
      glDeleteShader(shaders[i]);
    }
    glDeleteProgram(program_);
  }
  program_ = 0;
  source_.clear();
  log_.clear();
  uniforms_.clear();
}

bool Shader::compile(const string& defines,
                     const string& vertex_shader,
                     const string& fragment_shader) {
  clear();

  GLuint p = glCreateProgram();

  GLuint v = glCreateShader(GL_VERTEX_SHADER);
  {
    const char* srcs[2] = {defines.c_str(), vertex_shader.c_str()};
    glShaderSource(v, 2, srcs, 0);
    glCompileShader(v);
  }

  char log[2048];
  int logLength;

  glGetShaderInfoLog(v, sizeof(log), &logLength, log);
  if (logLength) {
    log_.append("--vertex:\n");
    log_.append(log);
  }

  GLuint f = glCreateShader(GL_FRAGMENT_SHADER);
  {
    const char* srcs[2] = {defines.c_str(), fragment_shader.c_str()};
    glShaderSource(f, 2, srcs, 0);
    glCompileShader(f);
  }

  glGetShaderInfoLog(f, sizeof(log), &logLength, log);
  if (logLength) {
    log_.append("--fragment:\n");
    log_.append(log);
  }

  glAttachShader(p, v);
  glAttachShader(p, f);
  glLinkProgram(p);

  glGetProgramInfoLog(p, sizeof(log), &logLength, log);
  if (logLength) {
    log_.append("--link:\n");
    log_.append(log);
  }

  {
	  // Capture active uniforms
	  GLint nUniforms = 0, maxLen = 0;
	  glGetProgramiv(p, GL_ACTIVE_UNIFORM_MAX_LENGTH, &maxLen);
	  glGetProgramiv(p, GL_ACTIVE_UNIFORMS, &nUniforms);

	  GLchar* name = (GLchar*)malloc(maxLen);

	  GLint size, location;
	  GLsizei written;
	  GLenum type;

	  printf(" Location | Type | Name\n");
	  printf("------------------------------------------------\n");
	  for( int i = 0; i < nUniforms; ++i ) {
		glGetActiveUniform(p, i, maxLen, &written,
				   &size, &type, name );
		location = glGetUniformLocation(p, name);
		printf(" %-8d | %-5x| %s\n", location, type, name);

		string sname(name);
		string stype;
		switch (type) {
			case 0x1404: stype.assign("int"); break;
			case 0x1406: stype.assign("float"); break;
			case 0x140a: stype.assign("double"); break;
			case 0x8b50: stype.assign("vec2"); break;
			case 0x8b51: stype.assign("vec3"); break;
			case 0x8b52: stype.assign("vec4"); break;
			case 0x8b53: stype.assign("ivec2"); break;
			case 0x8b54: stype.assign("ivec3"); break;
			case 0x8b55: stype.assign("ivec4"); break;
			case 0x8ffc: stype.assign("dvec2"); break;
			case 0x8ffd: stype.assign("dvec3"); break;
			case 0x8ffe: stype.assign("dvec4"); break;
			case 0x8b5e: stype.assign("sampler2d"); break;
		}
		if (!stype.empty() && sname.find("[") == string::npos) {
			uniforms_.append("uniform " + stype + " " + sname + ";\n");
		}
	  }

	  free(name);
  }

  program_ = p;

  source_.assign(defines);
  source_.append(fragment_shader);

  GLint status;
  glGetProgramiv(p, GL_LINK_STATUS, &status);

  return status == GL_TRUE;
}
