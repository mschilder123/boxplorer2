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

double fabs(double f) {
  return sign(f)*f;
}

dvec2 complexMul(dvec2 a, dvec2 b) {
  return dvec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

// MaddHattPatt youtube inspired

vec3 HSVtoRGB(vec3 hsv) {
  vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
  vec3 p = abs(fract(hsv.xxx + k.xyz) * 6.0 - k.www);
  return hsv.z * mix(k.xxx, clamp(p - k.xxx, 0.0, 1.0), hsv.y);
}

vec3 SimpleSmoothedLoopedHue(int iteration, double dist, dvec2 z, int maxIters) {
  float smoothIter = iteration + 1.0 - log2(log2(float(dot(z,z))));
  float t = float(smoothIter) / float(maxIters);

  vec3 color = vec3(0.0);

  if (iteration < maxIters) {
    vec3 hsv = vec3(t * 1.0 - (0.1 * time), 10.0 * dist, 1.0);
    color = HSVtoRGB(hsv);
  }

  return color;
}

#define TrapVector par[4]

vec3 getColor2D(dvec2 c, int maxIters) {
  dvec2 z = dvec2(0.0);
  int i = 0;
  double dist = 1000.0;
  for (i = 0; i < maxIters; i++) {
    z = complexMul(z,z) + c;
    dvec2 delta = z - TrapVector.xy;
    dist = min(dist, dot(delta, delta));
    double d = dot(z,z);
    if (d > 4.0) break;
  }
  return SimpleSmoothedLoopedHue(i, dist, z, maxIters);
}

void main() {
  dvec3 p = deye;

  // Intersect the view ray with the 2D plane at Z==0
  double totalD = -p.z / dir.z;
  p += totalD * dir;

  // Attempt at auto maxIters
  int maxIters = int(2.0*pow(log(float(3840.0 / totalD)), 1.8));

  vec3 col = vec3(0.0);
  if (totalD > 0.0) col = getColor2D(dvec2(p), maxIters + iters);

  // Write zBuffer and pixel
  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / float(clamp(totalD/length(dir), zNear, zFar)));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
