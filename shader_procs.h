#ifndef SHADER_PROCS_H
#define SHADER_PROCS_H

// Resolve OpenGL shader procs.
bool enableShaderProcs(void);

#define NO_SDL_GLEXT
#include <SDL_opengl.h>
#include <SDL.h>

typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

#if (defined __APPLE__)
  #include <OpenGL/glu.h>
  #include <OpenGL/glext.h>
#if !defined(GL_DECLARE_ONLY)
#ifndef PFNGLUNIFORM1DPROC
void (*glUniform1d)(int, double) = 0;
void (*glUniform3dv)(int, int, double*) = 0;
#endif
#endif  // GL_DECLARE_ONLY
 #elif (defined __WIN32__)
#elif (defined __GNUC__)
  #include <GL/glx.h>
  #include <GL/glxext.h>
#else
#error "Unknown target"
#endif

#if !(defined __APPLE__)
#include "GL/glext.h"

#if defined(GL_DECLARE_ONLY)
#define DECLARE_GL_PROC(type, name) extern type name
#else
#define DECLARE_GL_PROC(type, name) type name = 0
#endif

DECLARE_GL_PROC(PFNGLCREATEPROGRAMPROC, glCreateProgram);
DECLARE_GL_PROC(PFNGLDELETEPROGRAMPROC, glDeleteProgram);
DECLARE_GL_PROC(PFNGLISPROGRAMPROC, glIsProgram);
DECLARE_GL_PROC(PFNGLCREATESHADERPROC, glCreateShader);
DECLARE_GL_PROC(PFNGLDELETESHADERPROC, glDeleteShader);
DECLARE_GL_PROC(PFNGLGETSHADERIVPROC, glGetShaderiv);
DECLARE_GL_PROC(PFNGLSHADERSOURCEPROC, glShaderSource);
DECLARE_GL_PROC(PFNGLCOMPILESHADERPROC, glCompileShader);
DECLARE_GL_PROC(PFNGLATTACHSHADERPROC, glAttachShader);
DECLARE_GL_PROC(PFNGLDETACHSHADERPROC, glDetachShader);
DECLARE_GL_PROC(PFNGLGETATTACHEDSHADERSPROC, glGetAttachedShaders);
DECLARE_GL_PROC(PFNGLLINKPROGRAMPROC, glLinkProgram);
DECLARE_GL_PROC(PFNGLUSEPROGRAMPROC, glUseProgram);
DECLARE_GL_PROC(PFNGLGETSHADERIVPROC, glGetProgramiv);
DECLARE_GL_PROC(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog);
DECLARE_GL_PROC(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog);
DECLARE_GL_PROC(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation);
DECLARE_GL_PROC(PFNGLUNIFORM1FPROC, glUniform1f);
DECLARE_GL_PROC(PFNGLUNIFORM1IPROC, glUniform1i);
DECLARE_GL_PROC(PFNGLUNIFORM1IVPROC, glUniform1iv);
DECLARE_GL_PROC(PFNGLUNIFORM2FVPROC, glUniform2fv);
DECLARE_GL_PROC(PFNGLUNIFORM3FVPROC, glUniform3fv);
DECLARE_GL_PROC(PFNGLUNIFORM4FVPROC, glUniform4fv);
DECLARE_GL_PROC(PFNGLGETACTIVEUNIFORMARBPROC, glGetActiveUniform);

DECLARE_GL_PROC(PFNGLUNIFORM1DPROC, glUniform1d);
DECLARE_GL_PROC(PFNGLUNIFORM3DVPROC, glUniform3dv);

DECLARE_GL_PROC(PFNGLGENERATEMIPMAPPROC, glGenerateMipmap);
DECLARE_GL_PROC(PFNGLGENFRAMEBUFFERSPROC, glGenFramebuffers);
DECLARE_GL_PROC(PFNGLDELETEFRAMEBUFFERSPROC, glDeleteFramebuffers);
DECLARE_GL_PROC(PFNGLBINDFRAMEBUFFERPROC, glBindFramebuffer);
DECLARE_GL_PROC(PFNGLGENRENDERBUFFERSPROC, glGenRenderbuffers);
DECLARE_GL_PROC(PFNGLDELETERENDERBUFFERSPROC, glDeleteRenderbuffers);
DECLARE_GL_PROC(PFNGLBINDRENDERBUFFERPROC, glBindRenderbuffer);
DECLARE_GL_PROC(PFNGLRENDERBUFFERSTORAGEPROC, glRenderbufferStorage);
DECLARE_GL_PROC(PFNGLFRAMEBUFFERRENDERBUFFERPROC, glFramebufferRenderbuffer);
DECLARE_GL_PROC(PFNGLFRAMEBUFFERTEXTURE2DPROC, glFramebufferTexture2D);
DECLARE_GL_PROC(PFNGLCHECKFRAMEBUFFERSTATUSPROC, glCheckFramebufferStatus);

//DECLARE_GL_PROC(PFNGLFRAMEBUFFERTEXTURELAYERPROC, glFramebufferTextureLayer);
//DECLARE_GL_PROC(PFNGLTEXSTORAGE3DPROC, glTexStorage3D);

#if defined(__WIN32__)
DECLARE_GL_PROC(PFNGLACTIVETEXTUREPROC, glActiveTexture);
#endif

#undef DECLARE_GL_PROC

#endif  //!__APPLE__

#ifndef GL_DECLARE_ONLY
bool enableShaderProcs(void) {
#ifndef __APPLE__

#define IMPORT_GL_PROC(type, name) \
  do { if (!(name = (type) SDL_GL_GetProcAddress(#name))) { \
    fprintf(stderr, "failed to import function " #name "\n"); \
  } } while (0)

  IMPORT_GL_PROC(PFNGLCREATEPROGRAMPROC, glCreateProgram);
  IMPORT_GL_PROC(PFNGLCREATESHADERPROC, glCreateShader);
  IMPORT_GL_PROC(PFNGLSHADERSOURCEPROC, glShaderSource);
  IMPORT_GL_PROC(PFNGLCOMPILESHADERPROC, glCompileShader);
  IMPORT_GL_PROC(PFNGLATTACHSHADERPROC, glAttachShader);
  IMPORT_GL_PROC(PFNGLGETSHADERIVPROC, glGetShaderiv);
  IMPORT_GL_PROC(PFNGLLINKPROGRAMPROC, glLinkProgram);
  IMPORT_GL_PROC(PFNGLUSEPROGRAMPROC, glUseProgram);
  IMPORT_GL_PROC(PFNGLGETSHADERIVPROC, glGetProgramiv);
  IMPORT_GL_PROC(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog);
  IMPORT_GL_PROC(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog);
  IMPORT_GL_PROC(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation);
  IMPORT_GL_PROC(PFNGLUNIFORM1FPROC, glUniform1f);
  IMPORT_GL_PROC(PFNGLUNIFORM1IPROC, glUniform1i);
  IMPORT_GL_PROC(PFNGLUNIFORM1IVPROC, glUniform1iv);
  IMPORT_GL_PROC(PFNGLUNIFORM2FVPROC, glUniform2fv);
  IMPORT_GL_PROC(PFNGLUNIFORM3FVPROC, glUniform3fv);
  IMPORT_GL_PROC(PFNGLUNIFORM4FVPROC, glUniform4fv);
  IMPORT_GL_PROC(PFNGLGETATTACHEDSHADERSPROC, glGetAttachedShaders);
  IMPORT_GL_PROC(PFNGLDELETEPROGRAMPROC, glDeleteProgram);
  IMPORT_GL_PROC(PFNGLDELETESHADERPROC, glDeleteShader);
  IMPORT_GL_PROC(PFNGLDETACHSHADERPROC, glDetachShader);
  IMPORT_GL_PROC(PFNGLISPROGRAMPROC, glIsProgram);
  IMPORT_GL_PROC(PFNGLGETACTIVEUNIFORMARBPROC, glGetActiveUniform);

  IMPORT_GL_PROC(PFNGLUNIFORM1DPROC, glUniform1d);
  IMPORT_GL_PROC(PFNGLUNIFORM3DVPROC, glUniform3dv);

  IMPORT_GL_PROC(PFNGLGENERATEMIPMAPPROC, glGenerateMipmap);
  IMPORT_GL_PROC(PFNGLGENFRAMEBUFFERSPROC, glGenFramebuffers);
  IMPORT_GL_PROC(PFNGLDELETEFRAMEBUFFERSPROC, glDeleteFramebuffers);
  IMPORT_GL_PROC(PFNGLBINDFRAMEBUFFERPROC, glBindFramebuffer);
  IMPORT_GL_PROC(PFNGLGENRENDERBUFFERSPROC, glGenRenderbuffers);
  IMPORT_GL_PROC(PFNGLDELETERENDERBUFFERSPROC, glDeleteRenderbuffers);
  IMPORT_GL_PROC(PFNGLBINDRENDERBUFFERPROC, glBindRenderbuffer);
  IMPORT_GL_PROC(PFNGLRENDERBUFFERSTORAGEPROC, glRenderbufferStorage);
  IMPORT_GL_PROC(PFNGLFRAMEBUFFERRENDERBUFFERPROC, glFramebufferRenderbuffer);
  IMPORT_GL_PROC(PFNGLFRAMEBUFFERTEXTURE2DPROC, glFramebufferTexture2D);
  IMPORT_GL_PROC(PFNGLCHECKFRAMEBUFFERSTATUSPROC, glCheckFramebufferStatus);

  //IMPORT_GL_PROC(PFNGLFRAMEBUFFERTEXTURELAYERPROC, glFramebufferTextureLayer);
  //IMPORT_GL_PROC(PFNGLTEXSTORAGE3DPROC, glTexStorage3D);

#if defined(__WIN32__)
  IMPORT_GL_PROC(PFNGLACTIVETEXTUREPROC, glActiveTexture);
#endif

#undef IMPORT_GL_PROC

#endif  // __APPLE__

  return true;
}
#endif  // GL_DECLARE_ONLY

#endif  // SHADER_PROCS_H
