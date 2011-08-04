/*
 * Minimalist glsl funcs and defs needed to compile
 * compliant shader code as plain C++.
 */
#ifndef _F_INCLUDE_GLSL_H__
#define _F_INCLUDE_GLSL_H__

#define _FAKE_GLSL_

namespace GLSL {  // Wrap in namespace so we collide less w/ globals.

#ifndef PI
#define PI 3.14159265358979324
#endif

class vec3 {
public:
vec3() : x(0), y(0), z(0) {}
vec3(float k) : x(k), y(k), z(k) {}
vec3(float xx, float yy, float zz) { x = xx; y = yy; z = zz; }
vec3(const vec3& b) : x(b.x), y(b.y), z(b.z) {}
vec3& operator=(const vec3& b) { x = b.x; y = b.y; z = b.z; return *this; }
vec3& operator*=(const float k) { x*= k; y*= k; z*= k; return *this; }
vec3& operator/=(const float k) { x/= k; y/= k; z/= k; return *this; }
vec3& operator+=(const float k) { x+= k; y+= k; z+= k; return *this; }
vec3& operator+=(const vec3& b) { x+= b.x; y+= b.y; z+= b.z; return *this; }
vec3 operator*(const vec3& b) const { return vec3(x*b.x, y*b.y, z*b.z); }
vec3 operator*(const float k) const { return vec3(x*k, y*k, z*k); }
vec3 operator-(const vec3& b) const { return vec3(x-b.x, y-b.y, z-b.z); }
vec3 operator-() const { return vec3(-x,-y,-z); }
vec3 operator+(const vec3& b) const { return vec3(x+b.x, y+b.y, z+b.z); }
vec3 operator/(const float k) const { return vec3(x/k, y/k, z/k); }
vec3 cross(const vec3& b) const {
  return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
}
void print(const char* a) const { printf("%s(%f,%f,%f)\n",a,x,y,z); }
vec3 yxz() const { return vec3(y,x,z); }
vec3 xzy() const { return vec3(x,z,y); }
vec3 zyx() const { return vec3(z,y,x); }
vec3 zxy() const { return vec3(z,x,y); }
float x; float y; float z;
};

class vec2 {
public:
  vec2(float a, float b) : x(a), y(b) {}
  vec3 xxy() { return vec3(x,x,y); }
  vec3 xyx() { return vec3(x,y,x); }
  vec3 yxx() { return vec3(y,x,x); }
float x; float y;
};

class vec4 {
public:
vec4() : x(0), y(0), z(0), w(0) {}
vec4(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
vec4(const vec3& v3, float ww) { x = v3.x; y = v3.y; z = v3.z; w = ww; }
vec3 xyz() const { return vec3(x,y,z); }
vec4& operator=(const vec4& a) { x = a.x; y = a.y; z = a.z; w = a.w; return *this; }
vec4& operator/=(const float k) { x/= k; y/= k; z/= k; w/= k; return *this; }
vec4& operator*=(const float k) { x*= k; y*= k; z*= k; w*= k; return *this; }
vec4 operator*(const vec4& b) const { return vec4(x*b.x, y*b.y, z*b.z, w*b.w); }
vec4 operator+(const vec4& b) const { return vec4(x+b.x, y+b.y, z+b.z, w+b.w); }
vec4 operator/(const float k) const { return vec4(x/k, y/k, z/k, w/k); }
float dot(const vec4& b) const { return x*b.x + y*b.y + z*b.z + w*b.w; }
float x; float y; float z; float w;
};

float mod(float a, float b) { return a - b*floor(a/b); }
vec3 mod(const vec3& a, float b) { return vec3(mod(a.x, b), mod(a.y, b), mod(a.z, b)); }
float max(float a, float b) { return a>b?a:b; }
float min(float a, float b) { return a<b?a:b; }
float dot(const vec3& a, const vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
float length(const vec3& a) { return sqrt(dot(a, a)); }
float abs(float a) { return fabs(a); }
vec3 abs(const vec3& a) { return vec3(abs(a.x), abs(a.y), abs(a.z)); }
float clamp(float v, float l, float h) { if (v < l) return l; if (v > h) return h; return v; }
vec3 clamp(const vec3& v, float l, float h) {
  return vec3(clamp(v.x, l, h), clamp(v.y, l, h), clamp(v.z, l, h)); }
vec3 normalize(const vec3& a) { return vec3(a) / length(a); }
float sign(float a) {if (a<0) return -1; else return 1;}
vec3 mix(const vec3& a, const vec3&b, float r) {
  float ra = 1-r;
  float rb = r;
  return vec3(a.x*ra + b.x*rb, a.y*ra + b.y*rb, a.z*ra + b.z*rb);
}
float radians(float degrees) { return (degrees * PI) / 180.0; }

class mat4 {
 public:
  mat4(float a1, float a2, float a3, float a4,
       float b1, float b2, float b3, float b4,
       float c1, float c2, float c3, float c4,
       float d1, float d2, float d3, float d4) :
    r1(a1,a2,a3,a4),
    r2(b1,b2,b3,b4),
    r3(c1,c2,c3,c4),
    r4(d1,d2,d3,d4)
  {}

  vec4 operator*(const vec4& b) {
    return vec4(r1.dot(b), r2.dot(b), r3.dot(b), r4.dot(b));
  }

  vec4 r1;
  vec4 r2;
  vec4 r3;
  vec4 r4;
};

// Rewrite accessors to be method calls.
#define xyz xyz()
#define xzy xzy()
#define zyx zyx()
#define yzx yzx()
#define yxz yxz()
#define yxx yxx()
#define xyx xyx()
#define xxy xxy()

// inout (reference params) need a rewrite.
#define INOUT(a,b) a & b

// These attributes are meaningless for now in C++.
#define varying
#define uniform

}  // namespace GLSL

#endif
