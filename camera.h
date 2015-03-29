#ifndef _F_CAMERA_H_
#define _F_CAMERA_H_

#include <math.h>
#include <string>

#include "params.h"

class KeyFrame {
 public:
  double v[16];  // view matrix
  double q[4];   // quaterion orientation
  double x[4];   // r4 splineable q

  // Declare common parameters.
#define PROCESS(a,b,c,d) a b;
  PROCESS_COMMON_PARAMS
#undef PROCESS

  // Par[] parameter array.
  float par[NUMPARS][3];  // min(|this|, |glsl|) gets sent to shader.

  // Shader other uniforms.
  int iunis[NUMPARS];
  float funis[NUMPARS];

  int n_iunis;
  int n_funis;
  bool isKey_;  // Whether this frame is actually a defined KeyFrame.

  KeyFrame();

  double distanceTo(const KeyFrame& other) const;

  double* right() { return &v[0]; }
  double* up() { return &v[4]; }
  double* ahead() { return &v[8]; }
  double* pos() { return &v[12]; }

  void setKey(bool key) { isKey_ = key; }
  bool isKey() const { return isKey_; }

  // Orthogonalize v[]
  void orthogonalize();

  // Rotate q[] by `deg` degrees around a normalized axis
  // and set rotation part of v[] to q[].
  // Behaves like `glRotate` without normalizing the axis.
  void rotate(double deg, double x, double y, double z);

  // Move camera in a direction relative to the view direction.
  // Behaves like `glTranslate`.
  void move(double x, double y, double z) {
    for (int i=0; i<3; i++) {
      pos()[i] += right()[i]*x + up()[i]*y + ahead()[i]*z;
    }
    iBackbufferCount = 0;
  }

  // Move camera in the normalized absolute direction `dir` by `len` units.
  void moveAbsolute(double* dir, double len) {
    for (int i=0; i<3; i++) {
      pos()[i] += len * dir[i];
    }
  }

  // Map a uniform name to a address within this.
  // Returns NULL on fail.
  void* map_address(const std::string& type, const std::string& name, int n);
};

#endif
