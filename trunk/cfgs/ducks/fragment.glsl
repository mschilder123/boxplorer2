// 2D Ducks
// shader parts from Syntopia's frag-ducks > https://github.com/Syntopia/Fragmentarium
// Samuel Monnier's 'Ducks' Fractal
// see http://www.algorithmic-worlds.net/blog/blog.php?Post=20110227
// bermarte: thanx 2 Syntopia/Marius/Rrrola
// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform int iters;  // Number of fractal iterations. {min=10 max=1000}

#define surfaceColor par[2]
#define R par[2].x // {min=0. max=50. step=.01}
#define G par[2].y // {min=0. max=50. step=.01}
#define B par[2].z // {min=0. max=50. step=.01}
#define C par[3].x // {min=0. max=50. step=.01}
//#define Divider par[0].x  // {min=0. max=50. step=.01}
//#define Power par[0].y  // {min=0. max=6. step=.01}
//#define Radius par[0].z // {min=0. max=5. step=.01}

#define PrItrs par[1].x  // {min=0. max=50. step=.01}
#define JuliaX par[6].x// {min=-10. max=20. step=.001}
#define JuliaY par[6].y// {min=-10. max=20. step=.001}


vec2 zlog(vec2 a) {
	float b =  atan(a.y,a.x);
	if (b>0.0) b-=2.0*3.1415;
	return vec2(log(length(a)),b);
}

//REMOVE or re-use later
vec2 complexMul(vec2 a, vec2 b) {
  return vec2(a.x*b.x - a.y*b.y,a.x*b.y + a.y * b.x);
}
//

vec2 c2 = vec2(JuliaX,JuliaY);

// Mandelbrot for c.x,c.y
vec3 getColor2D(vec2 c) {
  bool Julia=true;

  vec2 z = Julia ?  c : vec2(1.0,0.0);  

  float mean = 0.0;
  
  int i = 0;
  float dist = 10000.0;

  for (i = 0; i < iters; i++) {

     z = zlog(vec2(z.x,abs(z.y)))+(Julia ? c2 : c);
     if (float(i)>PrItrs) mean+=length(z);
  }
    //dist = min(dist, abs(length(z)-Radius));
  
  mean/=float(iters)-PrItrs;
  // from Inigo Quilez's Shader Toy:
	float co =   1.0 - log2(.5*log2(mean/C));
	return vec3( .5+.5*cos(6.2831*co+R),.5+.5*cos(6.2831*co + G),.5+.5*cos(6.2831*co +B) );
}

void main() {
  vec3 p = eye;

  // Intersect the view ray with the 2D plane at Z==0
  float totalD = -p.z / dir.z;
  p += totalD * dir;

  vec3 col = vec3(0.0);
  if (totalD > 0.0) col = getColor2D(vec2(p));

  // Write zBuffer and pixel
  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
