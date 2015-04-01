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

#include "setup.inc"
#line 24

// Colors. Can be negative or >1 for interestiong effects.
#endif  // _FAKE_GLSL_

uniform bool julia;
#define JuliaVector par[1]

#define MB_SCALE par[0].y  // {min=-3 max=3 step=.001}
#define MB_MINRAD2 par[0].x  // {min=0 max=5 step=.001}

#define DIST_MULTIPLIER par[8].z  // {min=.01 max=2.0 step=.01}
#define MAX_DIST 10.0

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

vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4) * 4.0,
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
  scale = vec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(-MB_SCALE)) / minRad2;
  
  float s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
  //for (int i=0; i<int(float(iters) * detail); i++) s*= ds;
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

#define JL par[0].z // {min=0 max=1 step=1}
float de_mandelbox(vec3 pos) {
  vec4 p = vec4(pos,1.0);  // p.w is the distance estimate
  vec4 P0;
  if (julia) P0 = vec4(JuliaVector, JL); else P0 = p;
  int it = iters; //int(float(iters)*detail);
if (rotationAngle == 0.)
  for (int i=0; i<it; i++) {
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + P0;
    //if (r2 > 1000.0) break;
  }
else
  for (int i=0; i<it; i++) {
    p = vec4(rotationMatrix * p.xyz, p.w);
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + P0;
    //if (r2 > 1000.0) break;
  }
#if 1
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
            - absScalePowIters) * DIST_MULTIPLIER;
#else
  return length(p.xyz) / abs(p.w);
#endif
}
DECLARE_DE(de_mandelbox)

// Compute the color at `pos`.
vec3 c_mandelbox(vec3 pos) {
  vec3 p = pos;
  vec3 P0;
  if (julia) P0 = JuliaVector; else P0 = p;
  float trap = 1.0;
if (rotationAngle == 0.)
  for (int i=0; i<color_iters; i++) {
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    float r2 = dot(p, p);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale.xyz + P0;
    trap = min(trap, r2);
  }
else
  for (int i=0; i<color_iters; i++) {
    p = vec3(rotationMatrix * p.xyz);
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    float r2 = dot(p, p);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale.xyz + P0;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 2.0);
  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}
DECLARE_COLORING(c_mandelbox)

float normal_eps = 0.000001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec2 Eps = vec2(0, max(normal_eps, d_pos));
#if 0
  return normalize(vec3(
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)
  ));
#else
  // Eiffie's denormal
  vec3 dn = vec3(d(pos-Eps.yxx),d(pos-Eps.xyx),d(pos-Eps.xxy));
  vec3 dp = vec3(d(pos+Eps.yxx),d(pos+Eps.xyx),d(pos+Eps.xxy));
  return (dp-dn)/(length(dp-vec3(d_pos))+length(vec3(d_pos)-dn));
#endif
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

// from http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl
vec3 rgb2hsv(vec3 c) {
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main() {
  vec3 eye_in, dp; 

  if (!setup_ray(eye, dir, eye_in, dp)) return;
#if 0
{
  vec2 p = -1. + 2.*gl_FragCoord.xy / vec2(xres, yres);
  float lp = length(p);
  detail = 1.0 - .5 * smoothstep(.1, .8, lp);
}
#endif

  init();

  float m_zoom = zoom * 0.5 / (xres * pow(detail, 3.0));
  float noise = noise(gl_FragCoord.xy * time);

  vec3 p = eye_in;
  float D = d(p);
  float side = sign(D);
  float totalD = side * noise * D;   // Randomize first step.

  // Intersect the view ray with the Mandelbox using raymarching.
  float m_dist = m_zoom * totalD;

//#define OVERSTEP
#ifdef OVERSTEP
  float os = 0.;
#define sf par[0].z  // {min=0 max=2 step=.01}
#endif  // oVERSTEP
  int steps;
  for (steps=0; steps<max_steps; steps++) {
    D = (side * d(p + totalD * dp));
    if (D < m_dist) break;
#ifdef OVERSTEP
    if (D < os) {
      // might have overstepped
      totalD -= os;
      steps--;
      continue;
    } else {
      os = sf * D;
      totalD += D + os;
    }
#else  // OVERSTEP
    totalD += D;
#endif  // OVERSTEP

    if (totalD > MAX_DIST) break;
    m_dist =  m_zoom * totalD;
  }

#if 1
  // If we got a hit, find desired distance to surface.
  // Likely our hit was lot closer than m_dist; make it approx. m_dist.
  if (D < m_dist) {
    for (int i = 0; i < 5; ++i) {
      totalD += D - m_dist;
      m_dist =  m_zoom * totalD;
      D = d(p + totalD * dp);
    }
  }
#endif

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
    vec3 n = normal(p, D);
    col = c(p);
#if 1
    col = blinn_phong(n, -dp,
                      normalize(eye_in+vec3(0,1,0)+dp),
                      col);

    col = mix(aoColor, col, ambient_occlusion(p, n, abs(D), side, m_dist));

    //vec3 hsv = rgb2hsv(col);
    //hsv.z /= clamp(.005 * (m_dist / (m_zoom * abs(speed))), 1.0, 10.0);
    //col = hsv2rgb(hsv);
#else
    // light things up a bit
    // carry a light to the right and up a bit from eye, propotional
    // to speed (i.e. size of head / ipd).
    vec3 light_d = vec3(gl_ModelViewMatrix * vec4(abs(speed) * 14.0, abs(speed) * 5.0, 0, 0));
    vec3 light = eye + light_d;
    // trace back from p to light to determine shadow.
    vec3 tolight = normalize(light - p);
    float lightDist = length(light - p);
    float shadow = 1.0f;  // fully lit
    float lightRayLength = abs(speed);
    for (int i = 0; i < steps; ++i) {
      float d = d(p + lightRayLength * tolight);
      lightRayLength += d;
      if (lightRayLength > lightDist - 2*m_dist) {
        // reached light
        break;
      }

      shadow = min(shadow, 1.0 * d / lightRayLength);

      if (d < m_dist) {
        // hit obstacle; fully shaded
        shadow = 0.0;
        break;
      }
    }

    col = blinn_phong(n, -dp,
                      tolight,
                      //normalize(eye_in+vec3(0,1,0)+dp),
                      col);

    //col = mix(aoColor, col, ambient_occlusion(p, n, abs(D), side, m_dist));

    vec3 hsv = rgb2hsv(col);
    hsv.z *= shadow;
    col = hsv2rgb(hsv);
#endif

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col,
                backgroundColor,
                clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
  } else {
    // Record a miss as zFar; might be interpreted by effect shaders.
    totalD = 65535.0*abs(speed);
  }

  // Glow is based on the number of steps.
  col = mix(col,
            glowColor,
            (float(steps)+noise)/float(max_steps) * glow_strength);

  write_pixel(dir, totalD, sqrt(detail) * col);
}
