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
  vec2 operator-(const vec2& b) const;
  vec2 operator+(const vec2& b) const;
  vec2 operator*(const float b) const;
  vec2 operator*(const vec2& b) const;
  vec2 operator/(const vec2& b) const;
  vec2 operator+=(const vec2& b);
  vec2 operator-=(const vec2& b);
  float x,y;
};

class dvec3 {
 public:
  dvec3();
  dvec3(const double* v);
  dvec3(const dvec3& b);
  dvec3(double xx, double yy, double zz);
  dvec3 operator*(const double k) const;
  dvec3& operator=(const dvec3& b);
  dvec3& operator+=(const dvec3& b);
  dvec3 operator-(const dvec3& b) const;
  double x; double y; double z;
};

class vec3 {
 public:
  vec3();
  vec3(float k);
  vec3(const float* v);
  vec3(float xx, float yy, float zz);
  vec3(const vec2& b, float c);
  vec3(const vec3& b);
  vec3(const dvec3& b);
  vec3& operator=(const vec3& b);
  vec3& operator*=(const float k);
  vec3& operator/=(const float k);
  vec3& operator+=(const float k);
  vec3& operator+=(const vec3& b);
  vec3 operator*(const vec3& b) const;
  vec3 operator*(const float k) const;
  vec3 operator-(const vec3& b) const;
  vec3 operator-() const;
  vec3 operator+(const vec3& b) const;
  vec3 operator/(const float k) const;
  vec3 cross(const vec3& b) const;
  float dot(const vec3& b) const;
  double dot(const dvec3& b) const;
  void print(const char* a) const;
  vec3 yxz_() const;
  vec3 xzy_() const;
  vec3 zyx_() const;
  vec3 zxy_() const;
  vec2 xy_() const;
  vec2 xz_() const;
  float x; float y; float z;
};

class vec4 {
 public:
  vec4();
  vec4(float xx, float yy, float zz, float ww);
  vec4(const vec3& v3, float ww);
  vec4(const dvec3& v3, float ww);
  vec4(float v);
  vec3 xyz_() const;
  vec2 xy_() const;
  vec4& operator=(const vec4& a);
  vec4& operator/=(const float k);
  vec4& operator*=(const float k);
  vec4 operator*(const vec4& b) const;
  vec4 operator+(const vec4& b) const;
  vec4 operator/(const float k) const;
  float dot(const vec4& b) const;
  float x; float y; float z; float w;
};

class dvec4 {
 public:
  dvec4();
  dvec4(double xx, double yy, double zz, double ww);
  dvec4(const dvec3& v3, double ww);
  dvec4(double v);
  dvec3 xyz_() const;
  dvec4& operator=(const dvec4& a);
  dvec4& operator-=(const dvec4& a);
  dvec4& operator/=(const double k);
  dvec4& operator*=(const double k);
  dvec4 operator*(const dvec4& b) const;
  dvec4 operator+(const dvec4& b) const;
  dvec4 operator/(const double k) const;
  double dot(const dvec4& b) const;
  double x; double y; double z; double w;
};


vec2::vec2(float a) : x(a), y(a) {}
vec2::vec2(float a, float b) : x(a), y(b) {}
vec3 vec2::xxy_() const { return vec3(x,x,y); }
vec3 vec2::xyx_() const { return vec3(x,y,x); }
vec3 vec2::yxx_() const { return vec3(y,x,x); }
vec2 vec2::operator-(const vec2& b) const { return vec2(x-b.x, y-b.y); }
vec2 vec2::operator+(const vec2& b) const { return vec2(x+b.x, y+b.y); }
vec2 vec2::operator*(const float b) const { return vec2(x*b, y*b); }
vec2 vec2::operator*(const vec2& b) const { return vec2(x*b.x, y*b.y); }
vec2 vec2::operator/(const vec2& b) const { return vec2(x/b.x, y/b.y); }
vec2 vec2::operator-=(const vec2& b) { x-=b.x; y-=b.y; return *this; }
vec2 vec2::operator+=(const vec2& b) { x+=b.x; y+=b.y; return *this; }


vec3::vec3() : x(0), y(0), z(0) {}
vec3::vec3(const float* v) : x(v[0]), y(v[1]), z(v[2]) {}
vec3::vec3(float k) : x(k), y(k), z(k) {}
vec3::vec3(float xx, float yy, float zz) { x = xx; y = yy; z = zz; }
vec3::vec3(const vec2& b, float c) : x(b.x), y(b.y), z(c) {}
vec3::vec3(const vec3&b) : x(b.x), y(b.y), z(b.z) {}
vec3::vec3(const dvec3&b) : x(b.x), y(b.y), z(b.z) {}
vec3& vec3::operator=(const vec3& b) { x= b.x; y= b.y; z= b.z; return *this; }
vec3& vec3::operator*=(const float k) { x*= k; y*= k; z*= k; return *this; }
vec3& vec3::operator/=(const float k) { x/= k; y/= k; z/= k; return *this; }
vec3& vec3::operator+=(const float k) { x+= k; y+= k; z+= k; return *this; }
vec3& vec3::operator+=(const vec3& b) { x+= b.x; y+= b.y; z+= b.z; return *this; }
vec3 vec3::operator*(const vec3& b) const { return vec3(x*b.x, y*b.y, z*b.z); }
vec3 vec3::operator*(const float k) const { return vec3(x*k, y*k, z*k); }
vec3 vec3::operator-(const vec3& b) const { return vec3(x-b.x, y-b.y, z-b.z); }
vec3 vec3::operator-() const { return vec3(-x,-y,-z); }
vec3 vec3::operator+(const vec3& b) const { return vec3(x+b.x, y+b.y, z+b.z); }
vec3 vec3::operator/(const float k) const { return vec3(x/k, y/k, z/k); }
vec3 vec3::cross(const vec3& b) const { return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x); }
float vec3::dot(const vec3& b) const { return x*b.x + y*b.y + z*b.z; }
double vec3::dot(const dvec3& b) const { return ((double)x)*b.x + ((double)y)*b.y + ((double)z)*b.z; }
void vec3::print(const char* a) const { printf("%s(%f,%f,%f)\n",a,x,y,z); }
vec3 vec3::yxz_() const { return vec3(y,x,z); }
vec3 vec3::xzy_() const { return vec3(x,z,y); }
vec3 vec3::zyx_() const { return vec3(z,y,x); }
vec3 vec3::zxy_() const { return vec3(z,x,y); }
vec2 vec3::xy_() const { return vec2(x,y); }
vec2 vec3::xz_() const { return vec2(x,z); }

dvec3::dvec3(const dvec3&b) : x(b.x), y(b.y), z(b.z) {}
dvec3::dvec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
dvec3::dvec3(const double* v) : x(v[0]), y(v[1]), z(v[2]) {}
dvec3 dvec3::operator*(const double k) const { return dvec3(x*k, y*k, z*k); }
dvec3 dvec3::operator-(const dvec3& b) const { return dvec3(x-b.x, y-b.y, z-b.z); }
dvec3& dvec3::operator=(const dvec3& b) { x=b.x; y=b.y; z=b.z; return *this; }
dvec3& dvec3::operator+=(const dvec3& b) { x+=b.x;y+=b.y;z+=b.z;return *this;}


vec4::vec4() : x(0), y(0), z(0), w(0) {}
vec4::vec4(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
vec4::vec4(const vec3& v3, float ww) { x = v3.x; y = v3.y; z = v3.z; w = ww; }
vec4::vec4(float v) : x(v), y(v), z(v), w(v) {}
vec3 vec4::xyz_() const { return vec3(x,y,z); }
vec2 vec4::xy_() const { return vec2(x,y); }
vec4& vec4::operator=(const vec4& a) { x = a.x; y = a.y; z = a.z; w = a.w; return *this; }
vec4& vec4::operator/=(const float k) { x/= k; y/= k; z/= k; w/= k; return *this; }
vec4& vec4::operator*=(const float k) { x*= k; y*= k; z*= k; w*= k; return *this; }
vec4 vec4::operator*(const vec4& b) const { return vec4(x*b.x, y*b.y, z*b.z, w*b.w); }
vec4 vec4::operator+(const vec4& b) const { return vec4(x+b.x, y+b.y, z+b.z, w+b.w); }
vec4 vec4::operator/(const float k) const { return vec4(x/k, y/k, z/k, w/k); }
float vec4::dot(const vec4& b) const { return x*b.x + y*b.y + z*b.z + w*b.w; }


dvec4::dvec4() : x(0), y(0), z(0), w(0) {}
dvec4::dvec4(double xx, double yy, double zz, double ww) : x(xx), y(yy), z(zz), w(ww) {}
dvec4::dvec4(const dvec3& v3, double ww) { x = v3.x; y = v3.y; z = v3.z; w = ww; }
dvec4::dvec4(double v) : x(v), y(v), z(v), w(v) {}
dvec3 dvec4::xyz_() const { return dvec3(x,y,z); }
dvec4& dvec4::operator=(const dvec4& a) { x = a.x; y = a.y; z = a.z; w = a.w; return *this; }
dvec4& dvec4::operator-=(const dvec4& a) { x -= a.x; y -= a.y; z -= a.z; w -= a.w; return *this; }
dvec4& dvec4::operator/=(const double k) { x/= k; y/= k; z/= k; w/= k; return *this; }
dvec4& dvec4::operator*=(const double k) { x*= k; y*= k; z*= k; w*= k; return *this; }
dvec4 dvec4::operator*(const dvec4& b) const { return dvec4(x*b.x, y*b.y, z*b.z, w*b.w); }
dvec4 dvec4::operator+(const dvec4& b) const { return dvec4(x+b.x, y+b.y, z+b.z, w+b.w); }
dvec4 dvec4::operator/(const double k) const { return dvec4(x/k, y/k, z/k, w/k); }
double dvec4::dot(const dvec4& b) const { return x*b.x + y*b.y + z*b.z + w*b.w; }


float mod(float a, float b) { return a - b*floor(a/b); }
vec3 mod(const vec3& a, float b) { return vec3(mod(a.x, b), mod(a.y, b), mod(a.z, b)); }
//float max(float a, float b) { return a>b?a:b; }
double max(double a, double b) { return a>b?a:b; }
//float min(float a, float b) { return a<b?a:b; }
double min(double a, double b) { return a<b?a:b; }


float dot(const vec3& a, const vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double dot(const dvec3& a, const dvec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
float dot(const vec2& a, const vec2& b) { return a.x*b.x + a.y*b.y; }

vec3 cross(const vec3& a, const vec3& b) { return a.cross(b); }
float length(const vec3& a) { return sqrt(dot(a, a)); }
double length(const dvec3& a) { return sqrt(dot(a, a)); }
float length(const vec2& a) { return sqrt(dot(a, a)); }
float abs(float a) { return fabs(a); }
float fract(float a) { return fabs(a - floor(a)); }
double abs(double a) { return fabs(a); }
vec3 abs(const vec3& a) { return vec3(abs(a.x), abs(a.y), abs(a.z)); }
float clamp(float v, float l, float h) { if (v < l) return l; if (v > h) return h; return v; }
double clamp(double v, double l, double h) { if (v < l) return l; if (v > h) return h; return v; }
vec3 clamp(const vec3& v, float l, float h) {
  return vec3(clamp(v.x, l, h), clamp(v.y, l, h), clamp(v.z, l, h));
}
dvec3 clamp(const dvec3& v, double l, double h) {
  return dvec3(clamp(v.x, l, h), clamp(v.y, l, h), clamp(v.z, l, h));
}
vec3 clamp(const vec3& v, const vec3& l, const vec3& h) {
  return vec3(clamp(v.x, l.x, h.x), clamp(v.y, l.y, h.y), clamp(v.z, l.z, h.z));
}
vec2 clamp(const vec2& v, float l, float h) {
  return vec2(clamp(v.x, l, h), clamp(v.y, l, h));
}
vec3 normalize(const vec3& a) { return vec3(a) / length(a); }
float sign(float a) {if (a<0) return -1; else return 1;}
vec3 mix(const vec3& a, const vec3&b, float r) {
  float ra = 1-r;
  float rb = r;
  return vec3(a.x*ra + b.x*rb, a.y*ra + b.y*rb, a.z*ra + b.z*rb);
}
float radians(float degrees) { return float((degrees * PI) / 180.0); }
vec3 reflect(const vec3& d, const vec3& n) {
  return normalize(d - n*2*dot(n,d));
}
float sin(const float& a) { return (float)::sin(a); }
float cos(const float& a) { return (float)::cos(a); }
vec3 sin(const vec3& a) {
  return vec3(sin(a.x), sin(a.y), sin(a.z));
}
float floor(const float& a) { return (float)::floor(a); }

vec2 max(const vec2& a, const vec2& b) { return vec2(float(max(a.x, b.x)), float(max(a.y, b.y))); }
vec2 floor(const vec2& a) { return vec2(floor(a.x), floor(a.y)); }

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

class mat3 {
public:
  mat3(float a1, float a2, float a3,
       float b1, float b2, float b3,
       float c1, float c2, float c3) :
    r1(a1, a2, a3),
    r2(b1, b2, b3),
    r3(c1, c2, c3)
  {}

  dvec3 operator*(const dvec3& b) {
    return dvec3(r1.dot(b), r2.dot(b), r3.dot(b));
  }
  vec3 operator*(const vec3& b) {
    return vec3(r1.dot(b), r2.dot(b), r3.dot(b));
  }

  vec3 r1;
  vec3 r2;
  vec3 r3;
};

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
#define INOUT(a,b) a & b
#define IN_(a,b) a b
#define OUT_(a,b) a & b

// These attributes are meaningless for now in C++.
#define varying
#define uniform

}  // namespace GLSL

#endif
