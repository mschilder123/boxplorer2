#ifndef SHADER_PROCS_H
#define SHADER_PROCS_H

// Enable OpenGL 2.0 shader functions. Return 0 on error.
int enableShaderProcs(void);

////////////////////////////////

#define NO_SDL_GLEXT
#include <SDL/SDL_opengl.h>
#include <SDL/SDL.h>

#if (defined __APPLE__)
  #include <OpenGL/glu.h>
  #include <OpenGL/glext.h>
  int enableShaderProcs(void) { return 1; }
#elif (defined __WIN32__)
  #define GL_IMPORT_NEEDED
#elif (defined __linux__)
  #include <GL/glx.h>
  #include <GL/glxext.h>
  #define GL_IMPORT_NEEDED
#else
  int enableShaderProcs(void) { return 0; }
#endif


#ifdef GL_IMPORT_NEEDED

#include "GL/glext.h"

#define DECLARE_GL_PROC(type, name) \
  type name = 0
#define IMPORT_GL_PROC(type, name) \
  do { if (!(name = (type) SDL_GL_GetProcAddress(#name))) return 0; } while (0)

DECLARE_GL_PROC(PFNGLCREATEPROGRAMPROC, glCreateProgram);
DECLARE_GL_PROC(PFNGLDELETEPROGRAMPROC, glDeleteProgram);
DECLARE_GL_PROC(PFNGLISPROGRAMPROC, glIsProgram);
DECLARE_GL_PROC(PFNGLCREATESHADERPROC, glCreateShader);
DECLARE_GL_PROC(PFNGLDELETESHADERPROC, glDeleteShader);
DECLARE_GL_PROC(PFNGLSHADERSOURCEPROC, glShaderSource);
DECLARE_GL_PROC(PFNGLCOMPILESHADERPROC, glCompileShader);
DECLARE_GL_PROC(PFNGLATTACHSHADERPROC, glAttachShader);
DECLARE_GL_PROC(PFNGLDETACHSHADERPROC, glDetachShader);
DECLARE_GL_PROC(PFNGLGETATTACHEDSHADERSPROC, glGetAttachedShaders);
DECLARE_GL_PROC(PFNGLLINKPROGRAMPROC, glLinkProgram);
DECLARE_GL_PROC(PFNGLUSEPROGRAMPROC, glUseProgram);
DECLARE_GL_PROC(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog);
DECLARE_GL_PROC(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog);
DECLARE_GL_PROC(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation);
DECLARE_GL_PROC(PFNGLUNIFORM1FPROC, glUniform1f);
DECLARE_GL_PROC(PFNGLUNIFORM1IPROC, glUniform1i);
DECLARE_GL_PROC(PFNGLUNIFORM2FVPROC, glUniform2fv);
DECLARE_GL_PROC(PFNGLUNIFORM3FVPROC, glUniform3fv);

int enableShaderProcs(void) {
  IMPORT_GL_PROC(PFNGLCREATEPROGRAMPROC, glCreateProgram);
  IMPORT_GL_PROC(PFNGLCREATESHADERPROC, glCreateShader);
  IMPORT_GL_PROC(PFNGLSHADERSOURCEPROC, glShaderSource);
  IMPORT_GL_PROC(PFNGLCOMPILESHADERPROC, glCompileShader);
  IMPORT_GL_PROC(PFNGLATTACHSHADERPROC, glAttachShader);
  IMPORT_GL_PROC(PFNGLLINKPROGRAMPROC, glLinkProgram);
  IMPORT_GL_PROC(PFNGLUSEPROGRAMPROC, glUseProgram);
  IMPORT_GL_PROC(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog);
  IMPORT_GL_PROC(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog);
  IMPORT_GL_PROC(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation);
  IMPORT_GL_PROC(PFNGLUNIFORM1FPROC, glUniform1f);
  IMPORT_GL_PROC(PFNGLUNIFORM1IPROC, glUniform1i);
  IMPORT_GL_PROC(PFNGLUNIFORM2FVPROC, glUniform2fv);
  IMPORT_GL_PROC(PFNGLUNIFORM3FVPROC, glUniform3fv);
  IMPORT_GL_PROC(PFNGLGETATTACHEDSHADERSPROC, glGetAttachedShaders);
  IMPORT_GL_PROC(PFNGLDELETEPROGRAMPROC, glDeleteProgram);
  IMPORT_GL_PROC(PFNGLDELETESHADERPROC, glDeleteShader);
  IMPORT_GL_PROC(PFNGLDETACHSHADERPROC, glDetachShader);
  IMPORT_GL_PROC(PFNGLISPROGRAMPROC, glIsProgram);
  return 1;
}

#undef DECLARE_GL_PROC
#undef IMPORT_GL_PROC
#undef GL_IMPORT_NEEDED

#endif

#endif  // SHADER_PROCS_H
