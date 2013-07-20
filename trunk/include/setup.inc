uniform float focus;  // {min=-10 max=30 step=.01} Focal plane devation from 30x speed.

#ifndef OUT
#define OUT(a,b) out a b
#endif

bool setup_ray(vec3 eye, vec3 ray, OUT(vec3,eye_out), OUT(vec3,ray_out)) {
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
    //discard;  // Don't waste time on pixels we can't see.
    return false;
  }

  // Shift eye position, abs(speed) is half inter-occular distance.
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0.0, 0.0, 0.0));
  eye_out = eye + eye_d;

  // Note: no asymmetric frustum for Rift.
  ray_out = normalize(vec3(gl_ModelViewMatrix * vec4(p, 0.35, 0.0)));  // z value determines fov. Eyeballed.
#else
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed, 0, 0, 0));
#elif defined(ST_ANAGLYPH)
  float id = -1.0 + 2.0 * mod(gl_FragCoord.x + mod(gl_FragCoord.y, 2.0), 2.0);
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(id * speed, 0, 0, 0));
#else
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0, 0, 0));
#endif
  eye_out = eye + eye_d;
  // Construct asymmetric frustum.
  ray_out = normalize(ray * (focus + 30.0) * abs(speed) - eye_d);
#endif // ST_OCULUS
#else  // ST_NONE
  eye_out = eye;
  ray_out = normalize(ray);
#endif

  return true;
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
  gl_FragDepth = depth;
  gl_FragColor = vec4(finalCol, depth);
}