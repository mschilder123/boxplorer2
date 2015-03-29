// Super crude frame converter.

varying vec3 dir;

uniform float xres, yres, time;
uniform vec3 par[3];
#define WVector par[0]
#define AVector par[1]
#define BVector par[2]
uniform sampler2D iChannel0;
uniform int use_bg_texture;

vec4 color(vec2 position) {
  // position [0..1]
  vec2 p;
#if 0  // copy
  p = position;
#endif
#if 0  // o/u to sbs (l/r to l|r)
  if (position.x >= .5) {
    // right
    p = vec2((position.x - .5) * 2. + .0, (position.y) * .5 + .0);
  } else {
    // left
    p = vec2((position.x - .0) * 2. + .0, (position.y) * .5 + .5);
  }
#endif
#if 1  // r|l to cb (r|l to cb l|r)
  vec3 warp = par[0];
  vec2 asp = par[1].xy;
  vec2 asp2 = par[2].xy;
  if (position.x >= .5) {
    // right eye .75, .5
    p = position - vec2(.75, .5);
    p *= asp;
    float r2 = dot(p, p);
    float f = dot(warp, vec3(1.0, r2, r2*r2));
    p *= f;
    p *= asp2;
    p += vec2(.75, .5);
    if (p.x <= .5 || p.x >= 1. || p.y <= 0. || p.y >= 1.)
      return vec4(0.);
  } else {
    // left eye .25, .5
    p = position - vec2(.25, .5);
    p *= asp;
    float r2 = dot(p, p);
    float f = dot(warp, vec3(1.0, r2, r2*r2));
    p *= f;
    p *= asp2;
    p += vec2(.25, .5);
    if (p.x <= .0 || p.x >= .5 || p.y <= 0. || p.y >= 1.)
      return vec4(0.);
  }
  position = p;
#endif
  vec4 c = texture2D(iChannel0, p);
  return c;
}

void main() {
  vec2 pixelSize = vec2(1.0/xres, 1.0/yres);
  vec2 position = gl_FragCoord.xy * pixelSize;  // [0..1>
  vec4 c = color(position + pixelSize * .5);
  gl_FragColor = vec4(pow(c.rgb, vec3(1.0 / 2.2)), c.w);
}
