#pragma once

#include <string>

#ifndef ARRAYSIZE
#define ARRAYSIZE(x) (sizeof(x) / sizeof((x)[0]))
#endif

#define die(...) (fprintf(stderr, __VA_ARGS__), _exit(1), 1)

#define CHECK_STATUS(f, v)                                                     \
  {                                                                            \
    GLenum __s;                                                                \
    if ((__s = (f)) != (v)) {                                                  \
      printf("%s[%d] : %s() : %04x\n", __func__, __LINE__, #f, __s);           \
    }                                                                          \
  }

#define CHECK_ERROR CHECK_STATUS(glGetError(), GL_NO_ERROR)

#define CHECK_FRAMEBUFFER                                                      \
  CHECK_STATUS(glCheckFramebufferStatus(GL_FRAMEBUFFER),                       \
               GL_FRAMEBUFFER_COMPLETE)

// Global variables for paths
extern std::string BaseDir;
extern std::string WorkingDir;
extern std::string BaseFile;

// Utility functions
bool readFile(const std::string &name, std::string *content);
