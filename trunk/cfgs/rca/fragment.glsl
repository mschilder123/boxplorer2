// http://dmishin.blogspot.com/2013/11/the-single-rotation-rule-remarkably.html

varying vec3 dir;

uniform sampler2D iChannel0;
uniform float xres, yres, time;
uniform vec3 par[2];
uniform int use_bg_texture;
uniform int frameno;

#define decayR par[0].x  // {min=0 max=1 step=.0001}
#define decayG par[0].y  // {min=0 max=1 step=.0001}
#define decayB par[0].z  // {min=0 max=1 step=.0001}

#define random par[1].x  // {min=0 max=1 step=1}
#define direction par[1].y  // {min=-1 max=1 step=2}

#define backbuffer iChannel0
vec2 position;
vec2 pixelSize = vec2(1.0/xres, 1.0/yres);

float rand(vec2 co){
  return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void isAlive(float dx, float dy, inout int count, int factor) {
  float cur = texture2D(backbuffer, position + pixelSize*vec2( dx, dy )).r;
  // Aliveness is tracked in r channel at 1.0
  count += int(cur) * factor;
}

vec3 color(vec2 z) {
  // Ring o'fire; generate random live cells.
  if (use_bg_texture == 0 && random > 0.0) {
    if (length(z) < 0.1 && length(z) > 0.08) {
      // Within thin ring: random life!
      return (rand(time*z) < 0.5 ? vec3(1.0,0.0,0.0) : vec3(0.0));
    }
  }

  vec3 cur = texture2D(backbuffer, position).rgb;

  // Relative position of self in 2x2 cell.
  vec2 phase = -1.0 + 2.0 * floor(mod(gl_FragCoord.xy, 2.0));  // 1 or -1

  // Margolus toggle of which 2x2 to consider.
  phase *= -1.0 + 2.0 * floor(mod(float(frameno), 2.0));  // 1 or -1

  // Direction toggle.
  phase *= direction;

  int count = int(cur.r);  // self 1 or 0.

  // Count neighbours in relative 2x2.
  isAlive(phase.x, 0.0, count, 2);
  isAlive(phase.x, phase.y, count, 8);
  isAlive(0.0, phase.y, count, 4);
  
  // Rule: single rotate; get life if either count is 2 or 4,
  //       depending on rotation direction.
  // TODO: all Margolus rules?
  if (count == int(3.0 - direction * phase.x * phase.y)) {
    return vec3(1.0, 1.0, 1.0);  // rotate into new life
  } else if (count > 1) {  // more than 1 alive: stable
    return cur * vec3(1.0, par[0].y, par[0].z);
  } else {
    return cur * par[0];  // decay
  }
}

void main() {
  position = gl_FragCoord.xy * pixelSize;  // [0..1>
  vec2 z = -1.0 + 2.0 * position;          // <-1..1>
  z -= dir.xy;  // respond to mouse move, somewhat.
  gl_FragColor = vec4(color(z), 1.0);
}
