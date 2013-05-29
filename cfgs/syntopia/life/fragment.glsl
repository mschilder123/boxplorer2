// syntopia's Fragmentarium Life

varying vec3 dir;

uniform sampler2D bg_texture;
uniform float xres, yres, time;

#define backbuffer bg_texture
vec2 position;
vec2 pixelSize = vec2(1.0/xres, 1.0/yres);

float rand(vec2 co){
  return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void isAlive(float  dx, float  dy, inout int count) {
  vec4 v1 = texture2D( backbuffer,  position + pixelSize*vec2( dx, dy ) );
  if (v1.x==1.0) count++;
}

vec3 color(vec2 z) {
  // Ring o'fire
  if (length(z)<0.1 && length(z)>0.08) return (rand(time*z) < 0.5 ? vec3(1.0,0.0,0.0) : vec3(0.0));

  vec4 v1 = texture2D( backbuffer, position);
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
  
  // Rules
  if (alive==1) {
    if (neighbours<2) return vec3(v1.x*0.99,v1.y*0.99,v1.z);
    else if (neighbours<4) return vec3(1.0);
  } else {
    if (neighbours==3) return vec3(1.0);
  }

  return vec3(v1.x*0.95,v1.y*0.98,v1.z*0.999) ;
}

void main() {
  position = gl_FragCoord.xy * pixelSize;  // [0..1>
  vec2 z = -1.0 + 2.0 * position;          // <-1..1>
  z -= dir.xy;
  gl_FragColor = vec4(color(z), 1.0);
}
