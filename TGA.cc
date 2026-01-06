// Simple TGA handler. Only knows about 24bpp.

#include "TGA.h"

#include <SDL_opengl.h>
#include <stdio.h>

#include "glsl.h"

#if defined(_WIN32)
#pragma warning(disable : 4244)
#pragma warning(disable : 4996)
#endif

TGA::TGA() : data_(NULL) {}

TGA::~TGA() { delete[] data_; }

TGA::TGA(int width, int height)
    : width_(width), height_(height),
      data_(new unsigned char[width * height * 3]) {}

bool TGA::readFile(const char *filename) {
  FILE *f = NULL;
  bool result = false;
  while ((f = fopen(filename, "rb")) != NULL) {
    unsigned char header[18];
    if (fread(header, sizeof header, 1, f) != 1)
      break;
    unsigned char expected[18] = {0, 0, 2, 0, 0, 0, 0, 0,  0,
                                  0, 0, 0, 0, 0, 0, 0, 24, 0};
    expected[12] = header[12];
    expected[13] = header[13];
    expected[14] = header[14];
    expected[15] = header[15];
    if (memcmp(header, expected, sizeof header) != 0) {
      fprintf(stderr, "%s : unsupported TGA format, only 24bpp supported\n",
              __func__);
      break;
    }
    int width = header[13] * 256 + header[12];
    int height = header[15] * 256 + header[14];
    if (width > 32768 || height > 32768) {
      fprintf(stderr, "%s : oversized TGA image not supported\n", __func__);
      break;
    }
    unsigned char *data = new unsigned char[width * height * 3];
    if (fread(data, width * height * 3, 1, f) != 1) {
      fprintf(stderr, "%s : failed to load TGA pixel data\n", __func__);
      delete[] data;
      break;
    }
    delete[] data_;
    data_ = data;
    width_ = width;
    height_ = height;
    result = true;
    break;
  }
  if (f)
    fclose(f);
  return result;
}

bool TGA::writeFile(const char *filename) {
  const unsigned char header[18] = {0,
                                    0,
                                    2,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    (unsigned char)(width_ % 256),
                                    (unsigned char)(width_ / 256),
                                    (unsigned char)(height_ % 256),
                                    (unsigned char)(height_ / 256),
                                    24,
                                    0};
  FILE *f = NULL;
  bool result = false;
  if ((f = fopen(filename, "wb")) != NULL) {
    fwrite(header, 18, 1, f);
    fwrite(data_, 3, width_ * height_, f);
    result = true;
  }
  if (f)
    fclose(f);
  return result;
}

bool TGA::fromFramebuffer(int width, int height) {
  width_ = width;
  height_ = height;
  delete[] data_;
  data_ = new unsigned char[width_ * height_ * 3];
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  // glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, width_, height_, GL_BGR, GL_UNSIGNED_BYTE, data_);
  return true;
}

void TGA::set(int x, int y, GLSL::vec3 col) {
  data_[x * 3 + width_ * 3 * y + 0] = GLSL::clamp(256 * col.z, 0, 255); // B
  data_[x * 3 + width_ * 3 * y + 1] = GLSL::clamp(256 * col.y, 0, 255); // G
  data_[x * 3 + width_ * 3 * y + 2] = GLSL::clamp(256 * col.x, 0, 255); // R
}

bool TGA::get(int x, int y) {
  unsigned char r = data_[x * 3 + width_ * 3 * y + 2];
  return r == 255;
}

int TGA::width() { return width_; }
int TGA::height() { return height_; }
unsigned char *TGA::data() { return data_; }
