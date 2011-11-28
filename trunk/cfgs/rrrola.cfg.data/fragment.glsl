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
uniform float xres;

// Interactive parameters.
uniform vec3 par[10];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform float speed;

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

float de_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  vec4 scale = vec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  vec4 p = vec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  for (int i=0; i<iters; i++) {
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
           - pow(abs(MB_SCALE), float(1-iters))) * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

// Compute the color at `pos`.
vec3 c_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  vec3 scale = vec3(MB_SCALE, MB_SCALE, MB_SCALE) / minRad2;
  vec3 p = pos, p0 = p;
  float trap = 1.0;

  for (int i=0; i<color_iters; i++) {
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
  vec2 Eps = vec2(0, normal_eps);
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

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 20x speed.
void setup_stereo(inout vec3 eye_in, inout vec3 dp) {
#if !defined(ST_NONE)
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 4.0 * (fract(gl_FragCoord.y * 0.5) - .5) * abs(speed), 0, 0, 0));
#else
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0, 0, 0));
#endif
  eye_in = eye + eye_d;
  dp = normalize(dir * (focus + 20.0) * abs(speed) - eye_d);
#else  // ST_NONE
  eye_in = eye;
  dp = normalize(dir);
#endif
}

void main() {
  vec3 eye_in, dp; setup_stereo(eye_in, dp);

  float m_zoom = zoom * 0.5 / xres;

  float noise = pnoise(gl_FragCoord.xy);

  vec3 p = eye_in;
  float D = d(p);
  float side = sign(D);
  float totalD = side * D * noise;  // Randomize first step.

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

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
    vec3 n = normal(p, m_dist);
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

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
