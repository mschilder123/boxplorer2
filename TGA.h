// Simple TGA handler. Only knows about 24bpp.

#ifndef _F_INCLUDE_TGA_H__
#define _F_INCLUDE_TGA_H__

#include "glsl.h"
#include <string.h>

class TGA {
public:
  TGA();
  TGA(int width, int height);
  ~TGA();
  bool readFile(const char *filename);
  bool writeFile(const char *filename);
  bool fromFramebuffer(int width, int height);
  void set(int x, int y, GLSL::vec3 col);
  bool get(int x, int y);
  int width();
  int height();
  unsigned char *data();

private:
  int width_;
  int height_;
  unsigned char *data_;
};

#endif
