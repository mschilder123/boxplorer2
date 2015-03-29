// Lenord's Spudsville.
// Bits and pieces from knighty, rrrola, syntopia stitched together.

#define DIST_MULTIPLIER par[9].x  // {min=.01 max=1 step=.01}
#define MAX_DIST 10.0

#include "setup.inc"
#line 9

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;           // Distance at which raymarching stops.
uniform float ao_eps;             // Base distance at which ambient occlusion is estimated.
uniform float ao_strength;        // Strength of ambient occlusion.
uniform float glow_strength;      // How much glow is applied after max_steps.
uniform float dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters;                // Number of fractal iterations.
uniform int color_iters;          // Number of fractal iterations for coloring.
uniform int max_steps;            // Maximum raymarching steps.

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.1, 0.1, .4),
  aoColor = vec3(0, 0, 0);

#define Scale par[0].x  //{min=-5 max=4 step=.01}
#define fixedRadius2 par[1].x  //{min=.1 max=2.3 step=.01}
#define minRadius2 par[1].y  //{min=0 max=2.3 step=.01}
#define foldingLimit1 par[2].x  //{min=0 max=5 step=.01}
#define foldingLimit2 par[2].y  //{min=0 max=5 step=.01}
#define Power par[5].x  //{min=0.1 max=12.3 step=.1}
#define ZMUL par[5].y  //{min=-310 max=310 step=1}

void sphereFold(inout vec3 z, inout float dz) {
  float r2 = dot(z,z);
  if (r2 < minRadius2) {
    float temp = (fixedRadius2/minRadius2);
    z*= temp;
    dz*=temp;
  } else if (r2 < fixedRadius2) {
    float temp = (fixedRadius2/r2);
    z*= temp;
    dz*=temp;
  }
}

void boxFold1(inout vec3 z, inout float dz) {
  z = clamp(z, -foldingLimit1, foldingLimit1) * 2.0 - z;
}

void boxFold2(inout vec3 z, inout float dz) {
  z = clamp(z, -foldingLimit2, foldingLimit2) * 2.0 - z;
}

void powN2(inout vec3 z, float zr0, inout float dr) {
  float zo0 = asin( z.z/zr0 );
  float zi0 = atan( z.y,z.x );
  float zr = pow( zr0, Power-1.0 );
  float zo = zo0 * Power;
  float zi = zi0 * Power;
  dr = zr*dr*Power*length(vec3(1.0,1.0,ZMUL)/sqrt(3.0)) + 1.0;
  zr *= zr0;
  z = zr*vec3( cos(zo)*cos(zi), cos(zo)*sin(zi), ZMUL*sin(zo) );
}


float DE(vec3 z) {
  float dz = 1.0;
  float r = length(z);

  for (int n = 1; n < iters; ++n) {
    if (r >= 10.0) break;
    boxFold1(z,dz);
    sphereFold(z,dz);
    z *= Scale;
    dz *= abs(Scale);
    r = length(z);
  }

  if (r < 10.0) {
    boxFold2(z,dz);
    dz *= abs(Scale);
    r = length(z);
    powN2(z,r,dz);
    r = length(z);
  }

  return (r*log(r) / dz) * DIST_MULTIPLIER;
}

vec3 color(vec3 pos) {
  return vec3(0.7,0.6,0.4);
}

float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  //vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  vec2 Eps = vec2(0, max(d_pos, normal_eps));
  return normalize(vec3(
    -DE(pos-Eps.yxx)+DE(pos+Eps.yxx),
    -DE(pos-Eps.xyx)+DE(pos+Eps.xyx),
    -DE(pos-Eps.xxy)+DE(pos+Eps.xxy)
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
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side) {
  float ao_ed=DistAtp*ao_eps/min_dist;//Dividing by min_dist makes the AO effect independent from changing min_dist
  float ao = 1.0, w = ao_strength/ao_ed;//ps;
  float dist = 2.0 * ao_ed;//ps;

  for (int i=0; i<5; i++) {
    float D = side * DE(p + n*dist);
    ao -= (dist-D) * w;
    w *= 1.0/2.0;
    dist = dist*2.0 - ao_ed;//ps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

#define ULP 0.000000059604644775390625
float trace(inout vec3 p, in vec3 dp, inout float D, inout float totalD, in float side, in float MINDIST_MULT){
  // Intersect the view ray with the fractal using raymarching.
  // The distance field actually traced is the "calculated DE" minus (totalD * min_dist)
  // A perfect distance field have a gradient magnitude = 1. Assuming DE() gives a perfect DE, 
  //we have to multiply D with MINDIST_MULT in order to restore a gradient magnitude of 1
  int steps;
  for (steps=0; steps<max_steps &&
                abs(D)>max(totalD*8192.0*ULP,ULP) &&
                totalD < MAX_DIST; steps++) {
    totalD+=D;
    D = (side * DE(p + totalD * dp) - totalD * min_dist) * MINDIST_MULT;
  }
  p += totalD * dp;
  return float(steps);
}

float de_for_host(vec3 p) { return DE(p); }

void main() {
  vec3 eye_in, dp; 

  if (!setup_ray(eye, dir, eye_in, dp)) return;

  vec3 p = eye_in;

  float totalD = 0.0, D = DE(p);
  float side = sign(D);
  D = abs(D); 
  float MINDIST_MULT=1.0/(1.0+min_dist);
  D *= MINDIST_MULT;

  vec3 finalcol=vec3(0.);
  float refpart=1.0;
#define REFACTOR par[8].z  //{min=0 max=1 step=.01}
#define REFITER par[9].z
  bool cont=true;
  float firstD = 0.;  // save first step for depth buffer

  for (int ray = 0; ray < int(REFITER); ++ray) {
    if (!cont) break;

    float steps=trace(p, dp, D, totalD, side, MINDIST_MULT);
    if (ray == 0) firstD = totalD;
    vec3 col = backgroundColor;

    // We've got a hit or we're not sure.
    if (totalD < MAX_DIST) {
      float D1 = min_dist*.5*totalD;
      vec3 n = side * normal(p, max(256.0*ULP,D1));
      col = color(p);
      col = blinn_phong(n, -dp, normalize(eye_in+vec3(-1.0,-0.5,-0.7)), col);
      col = mix(aoColor, col, ambient_occlusion(p, n, D1, side));
  
      dp=reflect(dp,n);
      p -= totalD * dp;
      D = 9. * D1;

      if (D > max(totalD*8192.0*ULP,ULP)){
        col = mix(col, backgroundColor,
                  clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
      }
      col = mix(col, glowColor, steps/float(max_steps) * glow_strength);
    }else {
      cont=false;
    }

#if 1
    // Glow is based on the number of steps.
    col = mix(col, glowColor,
              (float(steps))/float(max_steps) * glow_strength);
#endif

    finalcol+=refpart*col;
    refpart*=REFACTOR;
  }

  write_pixel(dir, firstD, pow(finalcol, vec3(2.0)));
  //write_pixel(dir, firstD, finalcol);
}
