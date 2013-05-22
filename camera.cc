#include "camera.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>

#include "interpolate.h"

using namespace std;

KeyFrame::KeyFrame() {
  memset(v, 0, sizeof(v));
  memset(q, 0, sizeof(q));
  memset(x, 0, sizeof(x));
  // Initialize common parameters.
#define PROCESS(a,b,c,d) b = 0;
  PROCESS_COMMON_PARAMS
#undef PROCESS
  memset(par, 0, sizeof(par));
  isKey_ = false;
}

double KeyFrame::distanceTo(const KeyFrame& other) const {
  double delta[3] = { v[12]-other.v[12],
                      v[13]-other.v[13],
                      v[14]-other.v[14] };
  return sqrt(dot(delta, delta));
}

void KeyFrame::orthogonalize() {
  if (!normalize(ahead())) { ahead()[0]=ahead()[1]=0; ahead()[2]=1; }
  if (!normalize(up())) {  // Error? Make upDirection.z = 0.
    up()[2] = 0;
    if (fabs(ahead()[2]) == 1) {
      up()[0] = 0;
      up()[1] = 1;
    } else {
      up()[0] = -ahead()[1];
      up()[1] = ahead()[0];
      normalize(up());
    }
  }
  double l = dot(ahead(), up());
  for (int i=0; i<3; i++) up()[i] -= l*ahead()[i];
  // Compute rightDirection as a cross product of upDirection and direction.
  for (int i=0; i<3; i++) {
    int j = (i+1)%3, k = (i+2)%3;
    right()[i] = up()[j]*ahead()[k] - up()[k]*ahead()[j];
  }
  right()[3] = 0; up()[3] = 0; ahead()[3] = 0; pos()[3] = 1;
}

// Rotate the camera by `deg` degrees around a normalized axis.
// Behaves like `glRotate` without normalizing the axis.
void KeyFrame::rotate(double deg, double x, double y, double z) {
  quat2mat(this->q, this->v);
  double s = sin(deg*PI/180), c = cos(deg*PI/180), t = 1-c;
  double r[3][3] = {
    { x*x*t +   c, x*y*t + z*s, x*z*t - y*s },
    { y*x*t - z*s, y*y*t +   c, y*z*t + x*s },
    { z*x*t + y*s, z*y*t - x*s, z*z*t +   c }
  };
  for (int i=0; i<3; i++) {
    double c[3];
    for (int j=0; j<3; j++) c[j] = v[i+j*4];
    for (int j=0; j<3; j++) v[i+j*4] = dot(c, r[j]);
  }
  orthogonalize();
  mat2quat(this->v, this->q);
}

void* KeyFrame::map_address(const string& type, const string& name, int n) {
//  cerr << __FUNCTION__ << ": looking for " << type << " " << name << endl;

#define PROCESS(a,b,c,d) \
	if (name.compare(c) == 0) { \
      if (!d) return NULL; \
      if (type.compare(#a) == 0) return &this->b; \
      return NULL; \
    }
  PROCESS_COMMON_PARAMS
#undef PROCESS

  // TODO: for non-predefined uniforms, map them into couple of arrays,
  //       so we'd get automagic uniform discovery and linkage.

  return NULL;
}
