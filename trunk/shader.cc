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
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  GLuint f = glCreateShader(GL_FRAGMENT_SHADER);
  {
    const char* srcs[2] = {defines.c_str(), fragment_shader.c_str()};
    glShaderSource(f, 2, srcs, 0);
    glCompileShader(f);
  }

  glGetShaderInfoLog(f, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  glAttachShader(p, v);
  glAttachShader(p, f);
  glLinkProgram(p);

  glGetProgramInfoLog(p, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  {
	// Dump active uniforms
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
		printf(" %-8d | %-6x %s\n", location, type, name);
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
