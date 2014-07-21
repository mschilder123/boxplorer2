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
  vec4 cur = texture2D(backbuffer, position + pixelSize*vec2( dx, dy ));
  // Aliveness is tracked in r or a (except not a at first frame, to load img).
  float alive = max(cur.r, sign(float(frameno))*cur.a);
  count += int(alive) * factor;
}

vec4 color(vec2 z) {
  // Ring o'fire; generate random live cells.
  if (/*use_bg_texture == 0 && */random > 0.0) {
    if (length(z) < 0.02 /*&& length(z) > 0.08*/) {
      // Within tiny circle: random life!
      return (rand(time*z) < 0.5 ? vec4(1.0,0.0,0.0,1.0) : vec4(0.0));
    }
  }

  vec4 cur = texture2D(backbuffer, position);

  // Relative position of self in 2x2 cell.
  vec2 phase = -1.0 + 2.0 * floor(mod(gl_FragCoord.xy, 2.0));  // 1 or -1

  // Margolus toggle of which 2x2 to consider.
  phase *= -1.0 + 2.0 * floor(mod(float(frameno), 2.0));  // 1 or -1

  // Direction toggle.
  phase *= direction;

  float alive = max(cur.r, sign(float(frameno))*cur.a);
  int count = int(alive);

  // Count neighbours in relative 2x2.
  isAlive(phase.x, 0.0, count, 2);
  isAlive(0.0, phase.y, count, 4);
  isAlive(phase.x, phase.y, count, 8);
 
  // Bounce gas
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  //                  -  D  -  -  -  -  L  D  L  D  -  -  -  -  L  -
//const int rule[] = {0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1};
//int rule = 0xe968;
//int rotf = 0x0000;
//int rotb = 0x0000;

  // Billiard Ball Machine
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  //                  -  D  -  -  -  -  L  -  L  D  -  -  -  -  -  -
//const int rule[] = {0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1};
//int rule = 0xa9e8;
//int rotf = 0x0000;
//int rotb = 0x0000;

  // HPP Gas
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  //                  -  D  -  D  -  D  L  D  L  D  L  -  L  -  L  -
//const int rule[] = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1};
//int rule = 0xfd40;
//int rotf = 0x0000;
//int rotb = 0x0000;

  // SingleRotate
//const int rule[] = {0, 0, f, 1, b, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  //                  -  D  f  -  b  -  -  -  -  -  -  -  -  -  -  -
int rule = 0xaaa8;
int rotf = 0x0004;  // direction dependent toggle bits.
int rotb = 0x0010;

  // Scrambler
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
//const int rule[] = {0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1};
  //                  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
//int rule = 0xb315;
//int rotf = 0x0000;  // direction dependent toggle bits.
//int rotb = 0x0000;

#if 0
  if (fract(float(frameno) / 2.0) >= .5) {
    // Odd phase does Bounce gas.
    rule = 0xe968;
    rotf = 0;
    rotb = 0;
  }
#endif

  // add rotf if going forward.
  rule += rotf * int(clamp(direction * phase.x * phase.y, 0.0, 1.0));
  // add rotb if going backward.
  rule += rotb * int(clamp(-direction * phase.x * phase.y, 0.0, 1.0));

  if (fract(float(rule) / pow(2.0, float(count + 1))) < .5) {
    // dead: decay to black
    return cur * vec4(par[0], .5);
  } else if (fract(float(count) / 2.0) < .5) {
    const float t = .5 * .5;
    if (alive >= t) {
      // period 2: green
      return vec4(cur.r * .1, 1.0, cur.b * .45, 1.0);
    } else if (alive >= t * t) {
      // period 4: yellow
      return vec4(1.0, 1.0, cur.b * 1.0, 1.0);
    } else {
      // others: blue
      return vec4(cur.r * .3, cur.g * .4, 1.0, 1.0);
    }
  } else {
    // stable: climb to full red
    return vec4(clamp(cur.r + 0.001, 0.0, 1.0),
                cur.g * par[0].g, cur.b * par[0].b, 1.0);
  }
}

void main() {
  position = gl_FragCoord.xy * pixelSize;  // [0..1>
  vec2 z = -1.0 + 2.0 * position;          // <-1..1>
  z -= dir.xy;  // respond to mouse move, somewhat.
  gl_FragColor = color(z);
}
