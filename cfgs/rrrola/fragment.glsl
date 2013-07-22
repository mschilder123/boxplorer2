// Mandelbox shader by Rrrola
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal

#ifndef _FAKE_GLSL_

#define DECLARE_DE(x)
#define DECLARE_COLORING(x)
#define INOUT(a,b) inout a b

// distance estimator func
#ifndef d
#define d de_mandelbox // PKlein,combi,menger,mandelbox,ssponge
#endif

// surface coloring func
#ifndef c
#define c c_mandelbox  // PKlein,menger
#endif

#endif  // _FAKE_GLSL_

#define P0 p0                    // standard Mandelbox
//#define P0 vec4(par[1].x,par[1].y,par[2].y,1)  // Mandelbox Julia

#define MB_SCALE par[0].y  // {min=-3 max=3 step=.001}
#define MB_MINRAD2 par[0].x  // {min=0 max=5 step=.001}

#define DIST_MULTIPLIER par[8].z  // {min=.01 max=1.0 step=.01}
#define MAX_DIST 5.0

// Camera position and direction.
varying vec3 eye, dir;
varying float zoom;

uniform float xres, yres, speed;

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;           // Distance at which raymarching stops.
uniform float ao_eps;             // Base distance at which ambient occlusion is estimated.
uniform float ao_strength;        // Strength of ambient occlusion.
uniform float glow_strength;      // How much glow is applied after max_steps.
uniform float dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters;  // Number of fractal iterations. {min=0 max=500 step=1}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=0 max=500 step=1}
uniform int max_steps;  // Maximum raymarching steps. {min=0 max=1000 step=1}

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.03, 0.4, 0.4),
  aoColor = vec3(0, 0, 0);

float minRad2;
vec4 scale;
float absScalePowIters;
mat3 rotationMatrix;

void init() {
#define rotationVector par[3]
#define rotationAngle par[4].x  // { min=-5 max=5 step=.01}

  // compute couple of constants.
  minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  scale = vec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  
  float s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
  for (int i=0; i<iters; i++) s*= ds;
  absScalePowIters = s;
  
  float csat = cos(rotationAngle);
  float ssat = sin(rotationAngle);
  float usat = 1.0 - csat;
  vec3 u = normalize(rotationVector);
  rotationMatrix = mat3(
    u.x*u.x*usat + csat,     u.x*u.y*usat - u.z*ssat, u.x*u.z*usat + u.y*ssat,
    u.y*u.x*usat + u.z*ssat, u.y*u.y*usat + csat,     u.y*u.z*usat - u.x*ssat,
    u.z*u.x*usat - u.y*ssat, u.z*u.y*usat + u.x*ssat, u.z*u.z*usat + csat
    );
}

float de_mandelbox(vec3 pos) {
  vec4 p = vec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  for (int i=0; i<iters; i++) {
    p = vec4(rotationMatrix * p.xyz, p.w);
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
// if (r2 > 100.0) break;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
            - absScalePowIters) * 0.95 * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

// Compute the color at `pos`.
vec3 c_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  vec3 scale = vec3(MB_SCALE, MB_SCALE, MB_SCALE) / minRad2;
  vec3 p = pos, p0 = p;
  float trap = 1.0;
  for (int i=0; i<color_iters; i++) {
    p = vec3(rotationMatrix * p.xyz);
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    float r2 = dot(p, p);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);
  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}
DECLARE_COLORING(c_mandelbox)

float normal_eps = 0.000001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec2 Eps = vec2(0, max(normal_eps, d_pos));
  return normalize(vec3(
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)
  ));
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
}

// FAKE Ambient occlusion approximation.
// uses current distance estimate as first dist. the size of AO is independent from distance from eye
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side, float m_dist) {
  float ao_ed=ao_eps;
  float ao = 1.0, w = ao_strength/ao_ed;
  float dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-abs(D)) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

float hash( float n ) {
    return fract(sin(n)*5345.8621276);
}

float noise( in vec2 x ) {
    vec2 p = floor(x);
    vec2 f = fract(x);

    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*61.0;

    float res = mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                    mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y);

    return fract(res);
}

#include "setup.inc"

void main() {
  vec3 eye_in, dp; 

  if (!setup_ray(eye, dir, eye_in, dp)) {
    gl_FragColor = vec4(0.0);
    gl_FragDepth = 0.0;
    return;
  }

  init();

  float m_zoom = zoom * 0.5 / xres;
  float noise = noise(gl_FragCoord.xy / vec2(xres, yres));

  vec3 p = eye_in;
  float D = d(p);
  float side = sign(D);
  float totalD = side * D;   // Randomize first step.

  // Intersect the view ray with the Mandelbox using raymarching.
  float m_dist = m_zoom * totalD;
  int steps;
  for (steps=0; steps<max_steps; steps++) {
    D = (side * d(p + totalD * dp));
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist =  m_zoom * totalD;
  }

  // If we got a hit, find desired distance to surface.
  // Likely our hit was lot closer than m_dist; make it approx. m_dist.
  if (D < m_dist) {
    for (int i = 0; i < 5; ++i) {
      totalD += D - m_dist;
      m_dist =  m_zoom * totalD;
      D = d(p + totalD * dp);
    }
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
    vec3 n = normal(p, D/*m_dist*/);
    col = c(p);
    col = blinn_phong(n, -dp, normalize(eye_in+vec3(0,1,0)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, abs(D), side, m_dist));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor, clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);

  write_pixel(dir, totalD, col);
}
