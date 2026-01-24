// syntopia's Fragmentarium Life.
// coloring and tweaks by marius.

varying vec3 dir;

uniform sampler2D iBackbuffer;
uniform float xres, yres, time;
uniform vec3 par[2];

#define decayR par[0].x  // {min=0 max=1 step=.0001}
#define decayG par[0].y  // {min=0 max=1 step=.0001}
#define decayB par[0].z  // {min=0 max=1 step=.0001}
#define contagion par[1].x // {min=0 max=1 step=.0001}

vec2 position;
vec2 pixelSize = vec2(1.0/xres, 1.0/yres);

void Count(float dx, float dy, inout int count, inout int infected) {
  vec4 cur = texture2D(iBackbuffer, position + pixelSize*vec2(dx,dy));
  // Cell is alive if cur.a == 1.0
  count += int(cur.a);
  // Cell is infected if cur.g > .9 (border is solid 1.0)
  infected += int(cur.g + .1);
}

vec4 color() {
  vec4 cur = texture2D(iBackbuffer, position);
  int infected = int(cur.g + .1);
  int alive = int(cur.a);
  int neighbors = 0;
  
  // Count neighbors
  Count(1.0,0.0, neighbors, infected);
  Count(1.0,1.0, neighbors, infected);
  Count(1.0,-1.0, neighbors, infected);
  Count(0.0,1.0, neighbors, infected);
  Count(0.0,-1.0, neighbors, infected);
  Count(-1.0,1.0, neighbors, infected);
  Count(-1.0,0.0, neighbors, infected);
  Count(-1.0,-1.0, neighbors, infected);

  // Rules and coloring
  if (infected >= 1) {
    if (neighbors > 0) {
      // Spread infection, turn bright green infectuous
      return vec4(0,1,0,.5);
    } else {
      // Remain infectuous for a while
      return cur * vec4(decayR, contagion, decayB, .5);
    }
  }

  // Normal Conway rules
  if (alive == 1) {
    if (neighbors > 1 && neighbors < 4) {
      // Stay alive; turn red
      return vec4(clamp(cur.r + .0001, 0.0, 1.0),
                  cur.g * decayG,
                  cur.b * decayB,
                  1.0);
    }
  } else {
    if (neighbors == 3) {
      // Get born
      const float t = .5;
      if (cur.a >= t) {
        // Oscillate 1: mostly green
        return vec4(0.0,0.8,0.0,1.0);
      } else if (cur.a >= t * t) {
        // Oscillate 2: blue
        return vec4(0.0,0.0,1.0,1.0);
      } else {
        // Other: yellow-ish
        return vec4(0.9,0.8,0.0,1.0);
      }
    }
  }

  // Decaying
  return cur * vec4(par[0], .5);
}

void main() {
  position = gl_FragCoord.xy * pixelSize;  // [0..1>
  gl_FragColor = color();
}
