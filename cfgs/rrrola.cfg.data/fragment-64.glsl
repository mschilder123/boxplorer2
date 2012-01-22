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
uniform dvec3 deye;  // eye position in double precision

// Interactive parameters.
uniform vec3 par[10];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform double dspeed;

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

double de_mandelbox(dvec3 pos) {
  double minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  dvec4 scale = dvec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  dvec4 p = dvec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  double s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
  for (int i=0; i<iters; i++) {
    p = dvec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    double r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
	s *= ds;
  }
  return ((sqrt(dot(p.xyz, p.xyz)) - abs(MB_SCALE - 1.0)) / p.w
			- s) *.9;
         // - pow(abs(MB_SCALE), float(1-iters))) * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

// Compute the color at `pos`.
dvec3 c_mandelbox(dvec3 pos, double totalD) {
  double minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  dvec3 scale = vec3(MB_SCALE, MB_SCALE, MB_SCALE) / minRad2;
  dvec3 p = pos, p0 = p;
  double trap = 1.0;

  for (int i=0; i<color_iters; i++) {
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    double r2 = dot(p, p);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
//  dvec2 c = clamp(dvec2( 0.33*log(float(dot(p,p)))-1.0, sqrt(trap) ), 0.0, 1.0);
  dvec2 c = clamp(dvec2( (0.33*log(float(length(p)))-1.0)/totalD, sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
dvec3 normal(dvec3 pos, double d_pos) {
  dvec2 Eps = dvec2(0, d_pos);
  return normalize(dvec3(
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)
  ));
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
dvec3 blinn_phong(dvec3 normal, dvec3 view, dvec3 light, dvec3 diffuseColor) {
  dvec3 halfLV = normalize(light + view);
  double spe = max( dot(normal, halfLV), 0.0 );
  spe *= spe; // 2
  spe *= spe; // 4
  spe *= spe; // 8
  spe *= spe; // 16
  spe *= spe; // 32
  double dif = dot(normal, light) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
}

// FAKE Ambient occlusion approximation.
// uses current distance estimate as first dist. the size of AO is independent from distance from eye
double ambient_occlusion(dvec3 p, dvec3 n, double DistAtp, double side, double m_dist) {
  double ao_ed= m_dist; //ao_eps;
  double ao = 1.0, w = ao_strength/ao_ed;
  double dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    double D = side * d(p + n*dist);
    ao -= (dist-abs(D)) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 20x speed.
void setup_stereo(inout dvec3 eye_in, inout dvec3 dp) {
#if !defined(ST_NONE)
#if defined(ST_INTERLACED)
  dvec3 eye_d = dvec3(gl_ModelViewMatrix * dvec4( 4.0 * (fract(gl_FragCoord.y * 0.5) - .5) * abs(dspeed), 0, 0, 0));
#else
  dvec3 eye_d = dvec3(gl_ModelViewMatrix * dvec4(dspeed, 0, 0, 0));
#endif
  eye_in = deye + eye_d;
  dp = normalize(dir * (focus + 20.0) * abs(dspeed) - eye_d);
#else  // ST_NONE
  eye_in = deye;
  dp = normalize(dir);
#endif
}

void main() {
  dvec3 eye_in, dp; setup_stereo(eye_in, dp);

  double m_zoom = zoom * 0.5 / xres;

  double noise = pnoise(gl_FragCoord.xy);

  dvec3 p = eye_in;
  double D = d(p);
  double side = sign(D);
  double totalD = side * D * noise;  // Randomize first step.

  // Intersect the view ray with the Mandelbox using raymarching.
  double m_dist = m_zoom * totalD;
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
  dvec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
    dvec3 n = normal(p, m_dist);
    col = c(p, totalD);
    col = blinn_phong(n, -dp, normalize(eye_in+dvec3(0,1,0)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, abs(D), side, m_dist));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
 //   if (D > m_dist) {
 //     col = mix(col, backgroundColor, clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
 //   }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(float(col.x),float(col.y),float(col.z), depth);
}
