// syntopia's Fragmentarium Life.
// coloring and tweaks by marius.

varying vec3 dir;

uniform sampler2D iChannel0;
uniform float xres, yres, time;
uniform vec3 par[1];
uniform int use_bg_texture;

#define decayR par[0].x  // {min=0 max=1 step=.0001}
#define decayG par[0].y  // {min=0 max=1 step=.0001}
#define decayB par[0].z  // {min=0 max=1 step=.0001}

vec2 position;
vec2 pixelSize = vec2(1.0/xres, 1.0/yres);

float rand(vec2 co){
  return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void isAlive(float dx, float dy, inout int count) {
  vec4 cur = texture2D(iChannel0,  position + pixelSize*vec2( dx, dy ) );
  // Alive if cur.a == 1.0
  count += int(cur.a);
}

vec4 color(vec2 z) {
  // Ring o'fire
  if (use_bg_texture == 0) {
    if (length(z)<0.1 && length(z)>0.08) {
      return (rand(time*z) < 0.5 ? vec4(1.0,0.0,0.0,1.0) : vec4(0.0));
    }
  }

  vec4 cur = texture2D(iChannel0, position);
  int neighbours = 0;
  int alive = 0;

  isAlive(0.0,0.0,alive);
  
  // Count neighbours
  isAlive(1.0,0.0, neighbours);
  isAlive(1.0,1.0, neighbours);
  isAlive(1.0,-1.0, neighbours);
  isAlive(0.0,1.0, neighbours);
  isAlive(0.0,-1.0, neighbours);
  isAlive(-1.0,1.0, neighbours);
  isAlive(-1.0,0.0, neighbours);
  isAlive(-1.0,-1.0, neighbours);
  
  // Rules and coloring
  if (alive==1) {
    if (neighbours>1 && neighbours<4) {
      // Stay alive; turn red.
      return vec4(clamp(cur.r + .0001, 0.0, 1.0),
                  cur.g * par[0].g,
                  cur.b * par[0].b,
                  1.0);
    }
  } else {
    if (neighbours==3) {
      // born
      const float t = .5;
      if (cur.a >= t) {
        // Oscillate 1: green
        return vec4(0.0,1.0,0.0,1.0);
      } else if (cur.a >= t * t) {
        // Oscillate 2: blue
        return vec4(0.0,0.0,1.0,1.0);
      } else {
        // Other: yellow
        return vec4(1.0,1.0,0.0,1.0);
      }
    }
  }

  // Dead, decay
  return cur * vec4(par[0], .5);
}

void main() {
  position = gl_FragCoord.xy * pixelSize;  // [0..1>
  vec2 z = -1.0 + 2.0 * position;          // <-1..1>
  z -= dir.xy;
  gl_FragColor = color(z);
}
