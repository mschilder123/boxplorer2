#extension GL_ARB_gpu_shader_fp64 : enable

// 2D Mandelbrot
// shader parts from https://github.com/Syntopia/Fragmentarium

// Camera position and direction.
varying vec3 eye, dir;

#if 0
// wtf is wrong w/ AMD?
uniform dvec3 deye;  // eye position in double precision
#else
uniform double deyex, deyey, deyez;
dvec3 deye = dvec3(deyex, deyey, deyez);
#endif

// Interactive parameters.
uniform vec3 par[10];

uniform int iters;  // Number of fractal iterations. {min=10 max=10000}
uniform double dspeed;
uniform float time;

#define surfaceColor par[2]
#define surfaceColor2 par[3]
#define R par[2].x
#define G par[2].y
#define B par[2].z

double fabs(double f) {
  return sign(f)*f;
}

dvec2 complexMul(dvec2 a, dvec2 b) {
  return dvec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

#if 0
#define Divider par[0].x  // {min=0 max=50}
#define Power par[0].y  // {min=0 max=6 step=1e-4}
#define Radius par[0].z // {min=0 max=5 step=.01}

// Mandelbrot for c.x,c.y
vec3 getColor2D(dvec2 c, int maxIters) {
  dvec2 z = vec2(0.0,0.0);
  int i = 0;
  double dist = 10000.0;
  for (i = 0; i < maxIters; i++) {
    z = complexMul(z,z) + c;
    dist = min(dist, fabs(length(z)-Radius));
    if (dot(z,z)> 4.0) break;
  }
  if (i < maxIters) {
    // The color scheme here is based on one
    // from Inigo Quilez's Shader Toy:
    float co = float(i) + 1.0 - log2(.5*log2(float(dot(z,z))));
    co = sqrt(co/256.0);
    float co2 = dist * Divider;
    float fac = clamp(1.0/pow(co2,Power),0.0,1.0);
    return fac*vec3( .5+.5*cos(6.2831*co+R),
      .5+.5*cos(6.2831*co+G),
      .5+.5*cos(6.2831*co+B) );
  } else {
    return vec3(0.0);
  }
}
#endif

#if 1
// MaddHattPatt youtube inspired

vec3 HSVtoRGB(vec3 hsv) {
  vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
  vec3 p = abs(fract(hsv.xxx + k.xyz) * 6.0 - k.www);
  return hsv.z * mix(k.xxx, clamp(p - k.xxx, 0.0, 1.0), hsv.y);
}

vec3 SimpleSmoothedLoopedHue(int iteration, double dist, double z, int maxIters) {
  float smoothIter = iteration + 1.0 - log2(log2(float(length(z))));
  float t = float(smoothIter) / float(maxIters);

  vec3 color = vec3(0.0);

  if (iteration < maxIters) {
    vec3 hsv = vec3(t * 3.0 - (0.2 * time), 1.0 * dist, 1.0);
    color = HSVtoRGB(hsv);
  }
  return color;
}

#define TrapVector par[4]

vec3 getColor2D(dvec2 c, int maxIters) {
  dvec2 z = dvec2(0.0,0.0);
  int i = 0;
  double dist = 10000.0;
  for (i = 0; i < maxIters; i++) {
    z = complexMul(z,z) + c;
    dist = min(dist, length(z - TrapVector.xy));
    if (dot(z,z) > 16.0) break;
  }
  return SimpleSmoothedLoopedHue(i, dist, z, maxIters);
}

#endif


void main() {
  dvec3 p = deye;

  // Intersect the view ray with the 2D plane at Z==0
  double totalD = -p.z / dir.z;
  p += totalD * dir;

  // Attempt at auto maxIters
  int maxIters = int(2.0*pow(log(float(3840.0 / totalD)), 1.8));

  vec3 col = vec3(0.0);
  if (totalD > 0.0) col = getColor2D(dvec2(p), maxIters);

  // Write zBuffer and pixel
  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / float(clamp(totalD/length(dir), zNear, zFar)));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
