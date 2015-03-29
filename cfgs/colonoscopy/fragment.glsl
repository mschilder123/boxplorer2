// Mandelbox shader by Rrrola
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal

#define P0 p0                    // standard Mandelbox
//#define P0 vec4(par[1].x,par[1].y,par[2].y,1)  // Mandelbox Julia

#define MINRAD2 par[0].x  // {min=0 max=1 step=.001}
#define SCALE par[0].y  // {min=-3 max=3 step=.01}

#define fogStrength par[9].y  // {min=0 max=100 step=.25}
#define LoDbase par[8].x  // {min=1 max=20 step=.25}
#define LoDpow par[8].y  // {min=0 max=50 step=.25}

#define rotationAngle par[0].z  // {step=.001}
#define rotationVector par[2]

#define surfaceColor1 par[1]
#define surfaceColor2 par[3]
#define surfaceColor3 par[4]
#define glowColor par[5]
#define lightVector par[6]
#define backgroundColor par[7]

#define DIST_MULTIPLIER 0.5
#define MAX_DIST 4.0

// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;  // Distance at which raymarching stops. {min=1e-08 max=1e-04 step=1e-08}
uniform float ao_eps;  // Base distance at which ambient occlusion is estimated. {min=0 max=.001 step=.000001}
uniform float ao_strength;  // Strength of ambient occlusion. {min=0 max=.01 step=.0001}
uniform float glow_strength;  // How much glow is applied after max_steps. {min=0 max=10 step=.05}
uniform float dist_to_color;  // How is background mixed with the surface color after max_steps. {min=0 max=10 step=.01}

uniform float speed;  // {min=1e-06 max=.1 step=1e-06}

uniform int iters;  // Number of fractal iterations. {min=1 max=100}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=1 max=100}
uniform int max_steps;  // Maximum raymarching steps. {min=1 max=200}

uniform float xres, yres;

#include "setup.inc"

// Colors. Can be negative or >1 for interesting effects.
vec3 specularColor = vec3(.9, 0.8, 0.4);
vec3 aoColor = vec3(0, 0, 0);

// precomputed constants
float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
vec4 scale = vec4(SCALE, SCALE, SCALE, abs(SCALE)) / minRad2;
float absScalem1 = abs(SCALE - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(SCALE), float(1-iters));

// Rotate around vector
float csat = cos(rotationAngle);
float ssat = sin(rotationAngle);
float usat = 1.-cos(rotationAngle);
vec3 p0 = normalize(rotationVector);
mat3 RotationMatrix = mat3( p0.x*p0.x*usat + csat,      p0.x*p0.y*usat + p0.z*ssat, p0.x*p0.z*usat - p0.y*ssat,
                            p0.y*p0.x*usat - p0.z*ssat, p0.y*p0.y*usat + csat,      p0.y*p0.z*usat + p0.x*ssat,
                            p0.z*p0.x*usat + p0.y*ssat, p0.z*p0.y*usat - p0.x*ssat, p0.z*p0.z*usat + csat
                      );


// Compute the distance from `pos` to the Mandelbox.
float de_box(vec3 pos) {
  vec4 p = vec4(pos,1), p0 = p;  // p.w is the distance estimate

  for (int i=0; i<iters; i++) {
    // box folding: if (p>1) p = 2-p; else if (p<-1) p = -2-p;
//    p.xyz = abs(1.0+p.xyz) - p.xyz - abs(1.0-p.xyz);  // add;add;abs.add;abs.add (130.4%)
//    p.xyz = clamp(p.xyz*0.5+0.5, 0.0, 1.0) * 4.0 - 2.0 - p.xyz;  // mad.sat;mad;add (102.3%)
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;  // min;max;mad
//	 p.xyz *= RotationMatrix;
    // sphere folding: if (r2 < minRad2) p /= minRad2; else if (r2 < 1.0) p /= r2;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);  // dp3,div,max.sat,mul

    // scale, translate
    p = p*scale + P0;

   p.xyz *= RotationMatrix;
  }
  return ((length(p.xyz) - absScalem1) / p.w - AbsScaleRaisedTo1mIters) * DIST_MULTIPLIER;
}

// parameterized by iters
float de_box_iters(vec3 pos, int niters) {
  vec4 p = vec4(pos,1), p0 = p;  // p.w is the distance estimate

  for (int i=0; i<niters; i++) {
    // box folding: if (p>1) p = 2-p; else if (p<-1) p = -2-p;
//    p.xyz = abs(1.0+p.xyz) - p.xyz - abs(1.0-p.xyz);  // add;add;abs.add;abs.add (130.4%)
//    p.xyz = clamp(p.xyz*0.5+0.5, 0.0, 1.0) * 4.0 - 2.0 - p.xyz;  // mad.sat;mad;add (102.3%)
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;  // min;max;mad
//	p.xyz *= RotationMatrix;
    // sphere folding: if (r2 < minRad2) p /= minRad2; else if (r2 < 1.0) p /= r2;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);  // dp3,div,max.sat,mul
  
    // scale, translate
    p = p*scale + P0;

  p.xyz *= RotationMatrix;
  }
  return ((length(p.xyz) - absScalem1) / p.w - pow(abs(SCALE), float(1-niters))) * DIST_MULTIPLIER;
}

// Compute the color at `pos`.
vec3 color_box_iters(vec3 pos, int niters) {
  vec3 p = pos, p0 = p;
  float trap = 1.0;

  for (int i = 0; i < niters; i++) {
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
//	p.xyz *= RotationMatrix;
    p = p*scale.xyz + P0.xyz;
  p.xyz *= RotationMatrix;   // comment out for cool color interference, when rotated.
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}


//float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos, float normal_eps, int niters) {
  vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  return normalize(vec3(
  // 2-tap forward differences, error = O(eps)
//    -d_pos+de_box(pos+Eps.yxx),
//    -d_pos+de_box(pos+Eps.xyx),
//    -d_pos+de_box(pos+Eps.xxy)

  // 3-tap central differences, error = O(eps^2)
    -de_box_iters(pos-Eps.yxx, niters)+de_box_iters(pos+Eps.yxx, niters),
    -de_box_iters(pos-Eps.xyx, niters)+de_box_iters(pos+Eps.xyx, niters),
    -de_box_iters(pos-Eps.xxy, niters)+de_box_iters(pos+Eps.xxy, niters)

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
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor, float dist) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return clamp(dif*dist,0.,1.)*diffuseColor + (spe/dist)*specularColor;
}


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, int niters) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = de_box_iters(p + n*dist, iters);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_eps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

void main() {
  vec3 eye_in, dp; 
  if (!setup_ray(eye, dir, eye_in, dp)) return;

  vec3 p = eye_in;
  float totalD = 0.0, D = 0.0;

  int ncolor_iters = color_iters;

  // Intersect the view ray with the Mandelbox using raymarching.

  int steps;
  float m_dist = min_dist;

  for (steps = 0; steps < max_steps; steps++) {
    D = de_box(p + totalD * dp);
    if (D < m_dist) break;
    totalD += D;
    if (D > MAX_DIST) break;
    m_dist = min_dist * pow(totalD * LoDbase + 1.0, 100. * LoDpow);  // Vary LoD
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  float delta = m_dist / min_dist;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, D, m_dist * 10., iters);  // should look around in screen-space?

    col = color_box_iters(p, ncolor_iters);
    col = blinn_phong(n, -dp, normalize(eye_in+lightVector+dp), col, 1. /*clamp(delta/5., 1., 5.)*/);
    col = mix(aoColor, col, clamp(ambient_occlusion(p, n, iters), 0., 1.));

    // fog
    col = mix(col, backgroundColor, clamp(pow(1.+totalD*fogStrength, 2.)-1., 0., 1.));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor, clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  write_pixel(dir, totalD, col);
}
