// boxplorer2 inputs

uniform float aperture;  // {min=0 max=10 step=.01}

// Focal plane devation from 30x speed.
// Used for both asymmetric frustum and DoF,
// thus stereo convergence plane is in focus.
uniform float focus;  // {min=-20 max=30 step=.1}

// From vertex shader:
varying vec3 eye, dir;
varying float zoom;

// Uniforms:
uniform float xres, yres;

// Seconds since program / rendering sequence start.
uniform float time;

// Movement speed and half of virtual eye distance.
// Thus bigger virtual head means move faster.
uniform float speed;

#ifndef _FAKE_GLSL_
uniform sampler2D iBackbuffer;
uniform int iBackbufferCount;
uniform sampler2D iChannel0;
#endif  // _FAKE_GLSL_

uniform int frameno;

#if defined(ST_OCULUS)
// Physical ipd deviation from normal (i.e. me) for oculus user.
uniform float ipd;
#else
const
#endif
float detail = 1.0;  // full detail ratio.

#ifndef PI
const float PI = 3.14159265359;
#endif

#ifndef OUT_
#define OUT_(a,b) out a b
#endif

float de_for_host(vec3);
bool setup_ray(vec3 eye, vec3 ray, OUT_(vec3,eye_out), OUT_(vec3,ray_out)) {
#if defined(ST_DOME)
  vec2 a = -1. + 2. * gl_FragCoord.xy / vec2(xres, yres);  // [-1..1]
  if (length(a) > 1.) {
          gl_FragColor = vec4(0.);
          gl_FragDepth = 0.;
          return false;
  }
  a *= .5 * PI;  // [-PI/2..PI/2]
  float ca = cos(a.x);
  float sa = sin(a.x);
  mat3 ax = mat3(1., 0., 0.,
                 0., ca, sa,
                 0.,-sa, ca);
  float cb = cos(a.y);
  float sb = sin(a.y);
  mat3 ay = mat3(cb, 0.,-sb,
                 0., 1., 0.,
                 sb, 0., cb);
  float cc = cos(PI * .5);
  float sc = sin(PI * .5);
  mat3 az = mat3(cc, sc, 0.,
                 -sc, cc, 0.,
                 0., 0., 1.);
  vec3 x = vec3(0.,0.,1.) * ax * ay * az;
  eye_out = eye;
  ray_out = normalize(vec3(gl_ModelViewMatrix * vec4(x, 0.)));
  return true;
#endif      
#if defined(ST_SPHERICAL)
  vec2 a = -1. + 2. * gl_FragCoord.xy / vec2(xres, yres);
  a *= PI;
  a += PI;
  a.y += PI;
  a.y *= .5;
  float ca = cos(-a.x);
  float sa = sin(-a.x);
  mat3 ax = mat3(1., 0., 0.,
                 0., ca, sa,
                 0.,-sa, ca);
  float cb = cos(a.y);
  float sb = sin(a.y);
  mat3 ay = mat3(cb, 0.,-sb,
                 0., 1., 0.,
                 sb, 0., cb);
  float cc = cos(-PI * .5);
  float sc = sin(-PI * .5);
  mat3 az = mat3(cc, sc, 0.,
                 -sc, cc, 0.,
                 0., 0., 1.);
  vec3 x = vec3(0.,0.,1.) * ay * ax * az;
  eye_out = eye;
  ray_out = normalize(vec3(gl_ModelViewMatrix * vec4(x, 0.)));
  return true;
#endif
#if defined(ST_COMPUTE_DE_ONLY)
  // Caller is only interested in DE value as computed by shader.
  // Write it to output buffer, which better be FLOAT32.
  float dist = de_for_host(eye);
  gl_FragColor = vec4(dist);
  gl_FragDepth = dist;
  return false;
#else  // ST_COMPUTE_DE_ONLY
#if !defined(ST_NONE)
#if defined ST_OCULUS
  float halfx = xres / 2.0;

  vec2 q;
  if (sign(speed) < 0.0) {
    // left. 45 pixel shift towards center. Eyeballed.
    q = (gl_FragCoord.xy - vec2(-ipd * 10.0 + 20.0, 0.0)) / vec2(halfx, yres);
  } else {
    // right. 45 pixel shift towards center.
    q = (gl_FragCoord.xy - vec2(halfx + ipd * 10.0 - 20.0, 0.0)) / vec2(halfx, yres);
  }
  vec2 p = -1.0 + 2.0 * q;

  // Oculus barrel distort parameters.
  vec3 oculus_warp = vec3(1.0, 0.22, 0.24);  // k0, k1, k2
  vec2 oculus_scale = vec2(0.3, 0.35);  // x/y ratio eyeballed
  float r2 = dot(p, p);  // Radius squared, from center.
  float l1 = sqrt(r2);
  p *= oculus_scale * dot(oculus_warp, vec3(1.0, r2, r2*r2));
  float l2 = length(p);
  if (l2 > .5) { 
    // Don't waste time on pixels we can't see.
    gl_FragColor = vec4(0.5,0.,0.,0.);
    gl_FragDepth = 0.;
    return false;
  }
  if (l2 > .01) {
    // drop off detail towards periphery
    detail = 1. - .75 * smoothstep(.01, .5, l2);
  }

  // Shift eye position, abs(speed) is half inter-occular distance.
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0.0, 0.0, 0.0));
  eye_out = eye + eye_d;

  // Note: no asymmetric frustum for Rift.
  // z value (0.35) determines fov. Eyeballed.
  ray_out = normalize(vec3(gl_ModelViewMatrix * vec4(p, 0.35, 0.0)));
#else  // ST_OCULUS
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed, 0, 0, 0));
#elif defined(ST_ANAGLYPH)
  float id = -1.0 + 2.0 * mod(gl_FragCoord.x + mod(gl_FragCoord.y, 2.0), 2.0);
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(id * speed, 0, 0, 0));
#else  // ST_INTERLACED
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0, 0, 0));
#endif  // ST_INTERLACED
  eye_out = eye + eye_d;
  // Construct asymmetric frustum.
  ray_out = normalize(ray * (focus + 30.0) * abs(speed) - eye_d);
#endif // ST_OCULUS
#else  // ST_NONE
  eye_out = eye;
  ray_out = normalize(ray);
#endif

  return true;
#endif  // ST_COMPUTE_DE_ONLY
}

void write_pixel(vec3 ray, float z, vec3 finalCol) {
#if defined(ST_ANAGLYPH)
  float id = mod(gl_FragCoord.x + mod(gl_FragCoord.y, 2.0), 2.0);
  finalCol *= vec3(1.0 - id, id, id); 
#endif

  float zNear = abs(speed);

  // compute CoC, thin lens model
  float P = abs(focus + 30.0) * zNear;
  float D = max(zNear, z);;
  float A = aperture;
  float F = 8.*abs(speed); //~focalLength;
  float CoC = abs(A*(F*(P-D))/(D*(P-F)));

  float zFar = 65535.0 * zNear;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(D/length(ray), zNear, zFar));
  gl_FragColor = vec4(finalCol, clamp(CoC, float(0.0), float(1.0)));
  gl_FragDepth = depth;
}
