// Menger-sphere shader by marius. Parts by knighty, Rrrola.
// See http://www.fractalforums.com/3d-fractal-generation/revenge-of-the-half-eaten-menger-sponge/

#include "setup.inc"
#line 6

#define MAX_DIST 4.0

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;  // Distance at which raymarching stops. {min=1e-08 max=1e-04 step=1e-08}
uniform float ao_eps;  // Base distance at which ambient occlusion is estimated. {min=0 max=.001 step=.000001}
uniform float ao_strength;  // Strength of ambient occlusion. {min=0 max=.01 step=.0001}
uniform float glow_strength;  // How much glow is applied after max_steps. {min=0 max=10 step=.05}
uniform float dist_to_color;  // How is background mixed with the surface color after max_steps. {min=0 max=10 step=.05}

uniform int iters;  // Number of fractal iterations. {min=1 max=100}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=1 max=100}
uniform int max_steps;  // Maximum raymarching steps. {min=1 max=200}


// Colors. Can be negative or >1 for interesting effects.
#define backgroundColor par[7]
#define surfaceColor1 par[2]
#define surfaceColor2 par[3]
#define surfaceColor3 par[4]
#define specularColor par[5]
#define glowColor par[6]
#define lightVector par[1]
#define rotationVector par[9]

#define MOD par[8].x  // {min=1 max=5 step=.01}
#define SCALE par[8].y  // {min=1 max=5 step=.01}
#define INTRA_D par[8].z  // {min=.01 max=10.0 step=.01}
#define INITIAL_K par[0].z // {min=0 max=5 step=.01}
#define DE_EPS par[0].y  // {min=.0000001 max=.001 step=.0000001}

#define Angle par[0].z  // {min=-3 max=3 step=.001}
float csat = cos(Angle);
float ssat = sin(Angle);
float usat = 1.0-cos(Angle);
vec3 z0 = normalize(vec3(par[9].x, par[9].y, par[9].z));
mat3 RotationMatrix = mat3( z0.x*z0.x*usat + csat,      z0.x*z0.y*usat + z0.z*ssat, z0.x*z0.z*usat - z0.y*ssat,
                            z0.y*z0.x*usat - z0.z*ssat, z0.y*z0.y*usat + csat,      z0.y*z0.z*usat + z0.x*ssat,
			    z0.z*z0.x*usat + z0.y*ssat, z0.z*z0.y*usat - z0.x*ssat, z0.z*z0.z*usat + csat
		      );

vec3 aoColor = vec3(0, 0, 0);

//Functions to call.
#undef d
#define d d_SphereSponge
#define color color_SphereSponge

#if 1

float d_SphereSponge(vec3 pos) {
  float k = INITIAL_K;
  float d = -100.0;
  for(int i=0; i<iters; i++) {
    vec3 x = mod(pos * k, MOD) - .5 * MOD;
//    x *= RotationMatrix;
    float r = length(x);
    float d1 = (INTRA_D - r) / k;
    d = max(d, d1);
    k *= SCALE;
  }
  return d;
}

#else

// from http://fractal.io/
#define scale par[8].y  // {min=-10 max=10 step=.01}
#define offset par[8].x  // {min=-3 max=3 step=.01}
float d_SphereSponge(vec3 w) {
    w = (w * 0.5 + vec3(0.5)) * scale;  // scale [-1, 1] range to [0, 1]
    vec3 halfSpongeScale = vec3(.5) * scale;
    vec3 v = abs(w - halfSpongeScale) - halfSpongeScale;
    float d1 = max(v.x, max(v.y, v.z));     // distance to the box
    float d = d1;
    float p = 1.0;
    float md = 10000.0;
    vec3 cd = v;
    
    for (int i = 0; i < iters; i++) {
        vec3 a = mod(3.0 * w * p, 3.0);
        p *= 3.0;
        
        v = vec3(0.5) - abs(a - vec3(1.5)) + offset;
        v *= RotationMatrix;

        // distance inside the 3 axis aligned square tubes
        d1 = min(max(v.x, v.z), min(max(v.x, v.y), max(v.y, v.z))) / p;
        
        // intersection
        d = max(d, d1);
      
    }
    
    return d * 2.0 / scale;
}

#endif

vec3 color_SphereSponge(vec3 pos) {
#if 0
  float k = INITIAL_K;
  vec3 x;
  float trap = 1.0;
  for(int i = 0; i < color_iters; i++) {
    x = mod(pos * k, MOD) - .5 * MOD;
    float r = length(x);
    float d1 = (INTRA_D - r) / k;
    trap = min(trap, r);
    k *= SCALE;
  }
  return mix(surfaceColor1,
			 mix(surfaceColor2, surfaceColor3, length(x)),
			 trap);
#else
  return surfaceColor1;
#endif
}

float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  return normalize(vec3(
  // 2-tap forward differences, error = O(eps)
//    -d_pos+d(pos+Eps.yxx),
//    -d_pos+d(pos+Eps.xyx),
//    -d_pos+d(pos+Eps.xxy)

  // 3-tap central differences, error = O(eps^2)
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)

  // 4-tap forward differences, error = O(eps^3)
//    -2.0*d(pos-Eps.yxx)-3.0*d_pos+6.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
//    -2.0*d(pos-Eps.xyx)-3.0*d_pos+6.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
//    -2.0*d(pos-Eps.xxy)-3.0*d_pos+6.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)

  // 5-tap central differences, error = O(eps^4)
//    d(pos-Eps.zxx)-8.0*d(pos-Eps.yxx)+8.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
//    d(pos-Eps.xzx)-8.0*d(pos-Eps.xyx)+8.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
//    d(pos-Eps.xxz)-8.0*d(pos-Eps.xxy)+8.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)
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


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side) {
  float ao_ed=DistAtp*100.0;
  float ao = 1.0, w = ao_strength/ao_ed;//ps;
  float dist = 2.0 * ao_ed;//ps;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;//ps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

void main() {
  vec3 eye_in, dp;
  if (!setup_ray(eye, dir, eye_in, dp)) return;

  vec3 p = eye_in;

  float totalD = 0.0, D = 1000.0, extraD = 0.0, lastD;
  D = d(p + totalD * dp);
  float side = sign(D);
  D = D*side; D *= 0.95;
  totalD = D;

  // Intersect the view ray with the Mandelbox using raymarching.
  int steps;
  for (steps=0; steps<max_steps && D > max(DE_EPS*totalD,0.0000152587890625) && totalD < MAX_DIST; steps++) {
    D = side * d(p + totalD * dp);
    D = D - min_dist * (totalD + D) ;
    totalD += D;
  }

  p += totalD * dp;
  //normal_eps = max(normal_eps,1.0*D);
  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    D = max(D, max(DE_EPS*totalD,0.0000152587890625));
    vec3 n = side * normal(p, D);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(eye_in+lightVector+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, D, side));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > min_dist) {
      col = mix(col, backgroundColor, clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  // Write zBuffer and pixel
  write_pixel(dir, totalD, col);
}
