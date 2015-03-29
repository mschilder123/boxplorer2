#extension GL_ARB_gpu_shader_fp64 : enable

// Mandelbox shader by Rrrola
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal
// coloring and bits via fragmentarium
// double precision by marius

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
//#define P0 (p0-vec4(par[1],1.0))  // Mandelbox Julia
//#define P0 vec4(par[1],1.0)
//#define JuliaVector par[1]

#define MB_SCALE par[0].y  // {min=-3 max=3 step=.001}
#define MB_MINRAD2 par[0].x  // {min=0 max=5 step=.001}
#define detail par[0].z  // {min=.1 max=10 step=.01}

#define DIST_MULTIPLIER par[8].z  // {min=.01 max=1.0 step=.01}
#define MAX_DIST 15.0

// Camera position and direction.
varying vec3 eye, dir;
varying float zoom;
uniform float xres, yres;
uniform float aperture;  // {min=0.0 max=10.0 step=0.01}

#if 0
// wtf is wrong w/ AMD?
uniform dvec3 deye;  // eye position in double precision
#else
uniform double deyex, deyey, deyez;
dvec3 deye = dvec3(deyex, deyey, deyez);
#endif

// Interactive parameters.
uniform vec3 par[20];

uniform float min_dist;           // Distance at which raymarching stops.
uniform float ao_eps;             // Base distance at which ambient occlusion is estimated.
uniform float ao_strength;        // Strength of ambient occlusion.
uniform float glow_strength;      // How much glow is applied after max_steps.
uniform float dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform double dspeed;

uniform int iters;  // Number of fractal iterations. {min=0 max=500 step=1}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=0 max=500 step=1}
uniform int max_steps;  // Maximum raymarching steps. {min=0 max=1000 step=1}

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.06);

vec3 BaseColor = vec3(.3,.3,.3);
vec4 X = vec4(0.95, 0.64, 0.1, .9);
vec4 Y = vec4(0.89, 0.95, 0.75, .4);
vec4 Z = vec4(0.55, 0.16, 0.03, 1.0);
//vec4 X = vec4(0.95, 0.1, 0.1, .9);
//vec4 Y = vec4(0.1, 0.95, 0.1, .9);
//vec4 Z = vec4(0.1, 0.1, 0.95, .9);
vec4 R = vec4(0.12, 0.1, 0.2, 0.7);

vec3 specularColor = vec3(.9, 0.8, 0.4);
vec3 glowColor = vec3(0.3, 0.3, 0.3);
vec3 aoColor = vec3(0, 0, 0);

#define OrbitStrength par[9].x  //{min=0 max=1 step=.01}

double minRad2;
dvec4 scale;
double absScalePowIters;

mat3 rotationMatrix;
#define rotationVector par[3]
#define rotationAngle par[4].x  // { min=-5 max=5 step=.01}
#define rotationQuat par[3]  // that is par[3].xyz and par[4].x

void init() {
  // compute couple of constants.
  minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  scale = dvec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  
  double s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
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

double de_mandelbox(dvec3 pos) {
  dvec4 p = dvec4(pos, 1.0), p0 = p;  // p.w is the distance estimate
  for (int i=0; i<iters; i++) {
  	p.xyz = rotationMatrix * p.xyz;
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;

    double r2 = dot(p.xyz, p.xyz);

    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + P0;

	  if (r2 > 100.0) break;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
			- absScalePowIters) * 0.95 * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

// Compute the color at `pos`.
vec3 c_mandelbox(dvec3 pos, double totalD) {
  dvec4 p = dvec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  dvec4 orbitTrap = dvec4(10000.0);
  int i = 0;
  for (; i<color_iters; i++) {
    p.xyz = rotationMatrix * p.xyz;
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;

    double r2 = dot(p.xyz, p.xyz);
    orbitTrap = min(orbitTrap, abs(dvec4(p.xyz, r2)));

	  p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + P0;

	  if (r2 > 100.0) break;
  }

  orbitTrap.w = clamp(sqrt(orbitTrap.w)*float(color_iters)/float(i), 0.0, 1.0);
  vec3 orbitColor =
        X.xyz*X.w*float(orbitTrap.x) +
		Y.xyz*Y.w*float(orbitTrap.y) +
		Z.xyz*Z.w*float(orbitTrap.z) +
		R.xyz*R.w*float(orbitTrap.w)
		;
  return mix(BaseColor, 3.0*orbitColor,  OrbitStrength);
}


// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
dvec3 normal(dvec3 pos, double d_pos) {
  dvec2 Eps = dvec2(0, 3.0*d_pos);
  return normalize(
     dvec3(-abs(d(pos-Eps.yxx))+abs(d(pos+Eps.yxx)),
           -abs(d(pos-Eps.xyx))+abs(d(pos+Eps.xyx)),
           -abs(d(pos-Eps.xxy))+abs(d(pos+Eps.xxy))
		  )
     );
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(dvec3 normal, dvec3 view, dvec3 light, vec3 diffuseColor) {
  dvec3 halfLV = normalize(light + view);
  float spe = pow(float(max(dot(normal, halfLV), 0.0)), 16.0);
  float dif = float(dot(normal, light)) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
}

// FAKE Ambient occlusion approximation.
// uses current distance estimate as first dist. the size of AO is independent from distance from eye
float ambient_occlusion(dvec3 p, dvec3 n, double DistAtp) {
  double ao_ed = DistAtp;
  double ao = 1.0, w = ao_strength/ao_ed;
  double dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    double D = d(p + n*dist);
    ao -= (dist-abs(D)) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(float(ao), 0.0, 1.0);
}

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 30x speed.
bool setup_stereo(inout dvec3 eye_in, inout dvec3 dp) {
#if !defined(ST_NONE)
#if defined ST_OCULUS
  float halfx = xres / 2.0;

  vec2 q;
  if (sign(dspeed) < 0.0) {
    // left. 45 pixel shift towards center. Eyeballed.
    q = (gl_FragCoord.xy - vec2(focus + 45.0, 0.0)) / vec2(halfx, yres);
  } else {
    // right. 45 pixel shift towards center.
    q = (gl_FragCoord.xy - vec2(halfx - focus - 45.0, 0.0)) / vec2(halfx, yres);
  }
  vec2 p = -1.0 + 2.0 * q;

  // Oculus barrel distort parameters.
  vec3 oculus_warp = vec3(1.0, 0.22, 0.24);  // k0, k1, k2
  vec2 oculus_scale = vec2(0.3, 0.35);  // x/y ratio eyeballed
  float r2 = dot(p, p);  // Radius squared, from center.
  p *= oculus_scale * dot(oculus_warp, vec3(1.0, r2, r2*r2));
  if (dot(p, p) > 0.10) { 
    //discard;  // Don't waste time on pixels we can't see.
    return false;
  }

  // Shift eye position, abs(speed) is half inter-occular distance.
  dvec3 eye_d = dvec3(gl_ModelViewMatrix * dvec4(dspeed, 0.0, 0.0, 0.0));
  eye_in = deye + eye_d;

  // Note: no asymmetric frustum for Rift.
  dp = normalize(vec3(gl_ModelViewMatrix * vec4(p, 0.35, 0.0)));  // z value determines fov. Eyeballed.
#else
#if defined(ST_INTERLACED)
  dvec3 eye_d = dvec3(gl_ModelViewMatrix * dvec4( 4.0 * (fract(gl_FragCoord.y * 0.5) - .5) * dspeed, 0, 0, 0));
#else
  dvec3 eye_d = dvec3(gl_ModelViewMatrix * dvec4(dspeed, 0, 0, 0));
#endif
  eye_in = deye + eye_d;
  dp = normalize(dir * (focus + 30.0) * abs(dspeed) - eye_d);
#endif  // ST_OCULUS
#else  // ST_NONE
  eye_in = deye;
  dp = normalize(dir);
#endif
  return true;
}

void main() {
  dvec3 eye_in, dp;
  
  if (!setup_stereo(eye_in, dp)) {
    gl_FragColor = vec4(0);
	  gl_FragDepth = 0;
	  return;
  }

  init();

  double m_zoom = zoom * .5 /* sqrt(2)*/ / xres;

  double noise = pnoise(gl_FragCoord.xy);

  dvec3 p = eye_in;
  double D = d(p);
  double side = sign(D);
  double totalD = side * D * noise;  // Randomize first step.

  // Intersect the view ray with the Mandelbox using raymarching.
  double m_dist = m_zoom * totalD;
  int steps;
  for (steps=0; steps<max_steps; steps++) {
    D = (side * d(p + totalD * dp) - totalD * m_dist);
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist = m_zoom * totalD;
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
    col = c(p, totalD);
    dvec3 n = normal(p, m_dist/*abs(D)*/);
    col = blinn_phong(n, -dp, normalize(eye_in+dp+dvec3(0,-1,0)), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, 20.*m_dist));

#if 0
    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor, clamp(sqrt(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
#endif
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, (float(steps)+float(noise))/float(max_steps) * glow_strength);
  //col = mix(col, glowColor, clamp(float(totalD/dspeed/1000.0), 0.0, 1.0));

  double zNear = abs(dspeed);

  // compute CoC, thin lens model
  double P = abs(focus + 30.0) * zNear;
  D = totalD;
  double A = aperture;  //~aperture;
  double F = 8.*abs(dspeed); //~focalLength;
  double CoC = abs(A*(F*(P-D))/(D*(P-F)));

  double zFar = 65535.0 * zNear;
  double a = zFar / (zFar - zNear);
  double b = zFar * zNear / (zNear - zFar);
  float depth = float(a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragColor = vec4(col, float(clamp(CoC, 0.0, 1.0)));
  gl_FragDepth = depth;
}
