// Focal plane devation from 30x speed.
uniform float focus;  // {min=-10 max=30 step=.01}

#ifndef OUT_
#define OUT_(a,b) out a b
#endif

float de_for_host(vec3);
bool setup_ray(vec3 eye, vec3 ray, OUT_(vec3,eye_out), OUT_(vec3,ray_out)) {
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
    // Don't waste time on pixels we can't see.
    gl_FragColor = vec4(0.);
    gl_FragDepth = 0.;
    return false;
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
  float zFar = 65535.0 * zNear;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(z/length(ray), zNear, zFar));
  gl_FragColor = vec4(finalCol, depth);
  gl_FragDepth = depth;
}
