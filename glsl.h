/*
 * Minimalist glsl funcs and defs needed to compile
 * compliant shader code as plain C++.
 */
#pragma once
#define _FAKE_GLSL_

namespace GLSL { // Wrap in namespace so we collide less w/ globals.

#ifndef PI
#define PI 3.14159265358979324
#endif

class vec2;
class dvec2;
class vec3;
class dvec3;
class dvec4;

class vec2 {
public:
  vec2(float a);
  vec2(float a, float b);
  vec3 xxy_() const;
  vec3 xyx_() const;
  vec3 yxx_() const;
  vec2 operator-(const vec2 &b) const;
  vec2 operator+(const vec2 &b) const;
  vec2 operator*(const float b) const;
  vec2 operator*(const vec2 &b) const;
  vec2 operator/(const vec2 &b) const;
  vec2 operator+=(const vec2 &b);
  vec2 operator-=(const vec2 &b);

  float x;
  float y;
};

class dvec3 {
public:
  dvec3();
  dvec3(const double *v);
  dvec3(const dvec3 &b);
  dvec3(double xx, double yy, double zz);
  dvec3 operator*(const double k) const;
  dvec3 &operator=(const dvec3 &b);
  dvec3 &operator+=(const dvec3 &b);
  dvec3 operator-(const dvec3 &b) const;

  double x;
  double y;
  double z;
};

class vec3 {
public:
  vec3();
  vec3(float k);
  vec3(const float *v);
  vec3(float xx, float yy, float zz);
  vec3(const vec2 &b, float c);
  vec3(const vec3 &b);
  vec3(const dvec3 &b);
  vec3 &operator=(const vec3 &b);
  vec3 &operator*=(const float k);
  vec3 &operator/=(const float k);
  vec3 &operator+=(const float k);
  vec3 &operator+=(const vec3 &b);
  vec3 operator*(const vec3 &b) const;
  vec3 operator*(const float k) const;
  vec3 operator-(const vec3 &b) const;
  vec3 operator-() const;
  vec3 operator+(const vec3 &b) const;
  vec3 operator/(const float k) const;
  vec3 cross(const vec3 &b) const;
  float dot(const vec3 &b) const;
  double dot(const dvec3 &b) const;
  void print(const char *a) const;
  vec3 yxz_() const;
  vec3 xzy_() const;
  vec3 zyx_() const;
  vec3 zxy_() const;
  vec2 xy_() const;
  vec2 xz_() const;

  float x;
  float y;
  float z;
};

class vec4 {
public:
  vec4();
  vec4(float xx, float yy, float zz, float ww);
  vec4(const vec3 &v3, float ww);
  vec4(const dvec3 &v3, float ww);
  vec4(float v);
  vec3 xyz_() const;
  vec2 xy_() const;
  vec4 &operator=(const vec4 &a);
  vec4 &operator/=(const float k);
  vec4 &operator*=(const float k);
  vec4 operator*(const vec4 &b) const;
  vec4 operator+(const vec4 &b) const;
  vec4 operator/(const float k) const;
  float dot(const vec4 &b) const;

  float x;
  float y;
  float z;
  float w;
};

class dvec4 {
public:
  dvec4();
  dvec4(double xx, double yy, double zz, double ww);
  dvec4(const dvec3 &v3, double ww);
  dvec4(double v);
  dvec3 xyz_() const;
  dvec4 &operator=(const dvec4 &a);
  dvec4 &operator-=(const dvec4 &a);
  dvec4 &operator/=(const double k);
  dvec4 &operator*=(const double k);
  dvec4 operator*(const dvec4 &b) const;
  dvec4 operator+(const dvec4 &b) const;
  dvec4 operator/(const double k) const;
  double dot(const dvec4 &b) const;

  double x;
  double y;
  double z;
  double w;
};

class mat4 {
public:
  mat4(float a1, float a2, float a3, float a4, float b1, float b2, float b3,
       float b4, float c1, float c2, float c3, float c4, float d1, float d2,
       float d3, float d4);

  vec4 operator*(const vec4 &b);

  vec4 r1;
  vec4 r2;
  vec4 r3;
  vec4 r4;
};

class mat3 {
public:
  mat3(float a1, float a2, float a3, float b1, float b2, float b3, float c1,
       float c2, float c3);

  dvec3 operator*(const dvec3 &b);
  vec3 operator*(const vec3 &b);

  vec3 r1;
  vec3 r2;
  vec3 r3;
};

float mod(float a, float b);
vec3 mod(const vec3 &a, float b);
#undef max
// float max(float a, float b);
double max(double a, double b);
#undef min
// float min(float a, float b);
double min(double a, double b);

float dot(const vec3 &a, const vec3 &b);
double dot(const dvec3 &a, const dvec3 &b);
float dot(const vec2 &a, const vec2 &b);

vec3 cross(const vec3 &a, const vec3 &b);
float length(const vec3 &a);
double length(const dvec3 &a);
float length(const vec2 &a);
float abs(float a);
float fract(float a);
double abs(double a);
vec3 abs(const vec3 &a);
float clamp(float v, float l, float h);
double clamp(double v, double l, double h);
vec3 clamp(const vec3 &v, float l, float h);
dvec3 clamp(const dvec3 &v, double l, double h);
vec3 clamp(const vec3 &v, const vec3 &l, const vec3 &h);
vec2 clamp(const vec2 &v, float l, float h);
vec3 normalize(const vec3 &a);
float sign(float a);
vec3 mix(const vec3 &a, const vec3 &b, float r);
float radians(float degrees);
vec3 reflect(const vec3 &d, const vec3 &n);
float sin(const float &a);
float cos(const float &a);
vec3 sin(const vec3 &a);
float floor(const float &a);

vec2 max(const vec2 &a, const vec2 &b);
vec2 floor(const vec2 &a);

// Rewrite accessors to be method calls.
#define xyz xyz_()
#define xzy xzy_()
#define zyx zyx_()
#define yzx yzx_()
#define yxz yxz_()
#define yxx yxx_()
#define xyx xyx_()
#define xxy xxy_()
#define xy xy_()
#define xz xz_()

// inout (reference params) need a rewrite.
#define INOUT(a, b) a &b
#define IN_(a, b) a b
#define OUT_(a, b) a &b

// These attributes are meaningless for now in C++.
#define varying
#define uniform

} // namespace GLSL
