#pragma once

#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <string>

#ifndef ARRAYSIZE
#define ARRAYSIZE(x) (sizeof(x) / sizeof((x)[0]))
#endif

#define _pr_fmt(fmt) "%s:%d(%s) " fmt "\n", __FILE__, __LINE__, __func__
#define DIE(fmt, ...) (fprintf(stderr, _pr_fmt(fmt), ##__VA_ARGS__), exit(1), 1)

struct _LOG {
  _LOG(const char *file, int line, const char *func, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    printf("%s:%d(%s) ", file, line, func); // preample
    vprintf(fmt, args);
    va_end(args);
  }
  ~_LOG() { std::cout << std::endl; }

  template <typename T> friend _LOG &operator<<(_LOG &&ll, const T &t) {
    std::cout << t;
    return ll;
  }
};
#define DEBUG(fmt, ...) _LOG(__FILE__, __LINE__, __func__, fmt, __VA_ARGS__)

#define CHECK_STATUS(f, v)                                                     \
  {                                                                            \
    GLenum __s;                                                                \
    if ((__s = (f)) != (v)) {                                                  \
      DEBUG("%s() : %04x", #f, __s);                                           \
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
