// Simple TGA handler. Only knows about 24bpp.

#ifndef _F_INCLUDE_TGA_H__
#define _F_INCLUDE_TGA_H__

#include <string.h>

class TGA {
public:
  TGA() : data_(NULL) {}
  ~TGA() { delete[] data_; }
  TGA(size_t width, size_t height) :
      width_(width), height_(height),
      data_(new unsigned char[width*height*3]) {}
  bool readFile(const char* filename) {
    FILE* f = NULL;
    bool result = false;
    while ((f = fopen(filename, "rb")) != NULL) {
      unsigned char header[18];
      if (fread(header, sizeof header, 1, f) != 1) break;
      unsigned char expected[18] = {
        0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0
      };
      expected[12] = header[12]; expected[13] = header[13];
      expected[14] = header[14]; expected[15] = header[15];
      if (memcmp(header, expected, sizeof header) != 0) {
        fprintf(stderr, __FUNCTION__
                " : unsupported TGA format, only 24bpp supported\n");
        break;
      }
      size_t width = header[13] * 256 + header[12];
      size_t height = header[15] * 256 + header[14];
      if (width > 32768 || height > 32768) {
        fprintf(stderr, __FUNCTION__
                " : oversized TGA image not supported\n");
        break;
      }
      unsigned char* data = new unsigned char[width * height * 3];
      if (fread(data, width * height * 3, 1, f) != 1) {
        fprintf(stderr, __FUNCTION__
                " : failed to load TGA pixel data\n");
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
    if (f) fclose(f);
    return result;
  }
  bool writeFile(const char* filename) {
    const unsigned char header[18] = {
      0,0,2,0,0,0,0,0,0,0,0,0,
      (unsigned char)(width_%256),
      (unsigned char)(width_/256),
      (unsigned char)(height_%256),
      (unsigned char)(height_/256),24,0
    };
    FILE* f = NULL;
    bool result = false;
    if ((f = fopen(filename, "wb")) != NULL) {
      fwrite(header, 18, 1, f);
      fwrite(data_, 3, width_*height_, f);
      result = true;
    }
    if (f) fclose(f);
    return result;
  }
#ifdef GL_BGR
  bool readFramebuffer(size_t width, size_t height, int viewportOffset[2]) {
    width_ = width;
    height_ = height;
    delete[] data_;
    data_ = new unsigned char[width_ * height_ * 3];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    //glReadBuffer(GL_FRONT);
    glReadPixels(0,0,//viewportOffset[0], viewportOffset[1],
                 width_, height_,
                 GL_BGR, GL_UNSIGNED_BYTE,
                 data_);
    return true;
  }
#endif
#ifdef _FAKE_GLSL_
  void set(int x, int y, GLSL::vec3 col) {
    data_[x*3+width_*3*y+0] = GLSL::clamp(256*col.z, 0, 255);  // B
    data_[x*3+width_*3*y+1] = GLSL::clamp(256*col.y, 0, 255);  // G
    data_[x*3+width_*3*y+2] = GLSL::clamp(256*col.x, 0, 255);  // R
  }
#endif
  size_t width() { return width_; }
  size_t height() { return height_; }
  unsigned char* data() { return data_; }

private:
  size_t width_;
  size_t height_;
  unsigned char* data_;
};

#endif
