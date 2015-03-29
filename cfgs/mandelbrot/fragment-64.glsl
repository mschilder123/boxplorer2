#extension GL_ARB_gpu_shader_fp64 : enable
// 2D Mandelbrot
// shader parts from https://github.com/Syntopia/Fragmentarium

// Camera position and direction.
varying vec3 eye, dir;
#if 0
// wtf is wrong w/ AMD?
uniform dvec3 deye;  // eye position in double precision
#else
uniform double deyex, deyey, deyez;
dvec3 deye = dvec3(deyex, deyey, deyez);
#endif

// Interactive parameters.
uniform vec3 par[10];

uniform int iters;  // Number of fractal iterations. {min=10 max=10000}
uniform double dspeed;
uniform float time;

#define surfaceColor par[2]
#define surfaceColor2 par[3]
//#define R par[2].x
//#define G par[2].y
//#define B par[2].z

dvec2 complexMul(dvec2 a, dvec2 b) {
  return dvec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y * b.x);
}

#if 0
#define Divider par[0].x  // {min=0 max=50}
#define Power par[0].y  // {min=0 max=6 step=1e-4}
#define Radius par[0].z // {min=0 max=5}

// Mandelbrot for c.x,c.y
vec3 getColor2D(dvec2 c) {
  dvec2 z = vec2(0.0,0.0);
  int i = 0;
  double dist = 10000.0;
  for (i = 0; i < iters; i++) {
    z = complexMul(z,z) + c;
    if (dot(z,z)> 4.0) break;
    dist = min(dist, abs(length(z)-Radius));
  }
  if (i < iters) {
    // The color scheme here is based on one
    // from Inigo Quilez's Shader Toy:
    float co = float(i) + 1.0 - log2(.5*log2(float(dot(z,z))));
    co = sqrt(co/256.0);
    float co2 = dist * Divider;
    float fac = clamp(1.0/pow(co2,Power),0.0,1.0);
    return fac*vec3( .5+.5*cos(6.2831*co+R),
      .5+.5*cos(6.2831*co+G),
      .5+.5*cos(6.2831*co+B) );
  } else {
    return vec3(0.0);
  }
}
#endif

#if 0
// kali's grid deform coloring
// http://www.fractalforums.com/mandelbrot-and-julia-set/how-mandelbrot-deforms-a-grid/

#define Power par[0].y  // {min=0 max=6 step=1e-4}

vec3 getColor2D(dvec2 c) {
  dvec2 z = dvec2(0.0,0.0);
  for (int i = 0; i < iters; i++) {
    z = complexMul(z,z) + c;
    if (dot(z,z) > 4.0) break;
  }
  z.xy *= Power;
  float dx = fract(z.x);
  float dy = fract(z.y);
  if ((dx > .5 && dy > .5) ||
      (dx < .5 && dy < .5))
      return surfaceColor;
  else
      return surfaceColor2;
}
#endif

#if 1
// https://www.shadertoy.com/view/ldf3DN
vec3 getColor2D(dvec2 cc) {
  dvec2 z  = dvec2(0.0);
  dvec2 dz = dvec2(0.0);
  double trap1 = 0.0;
  double trap2 = 1e20;
  double co2 = 0.0;

#define trap2Vector par[4]
  dvec2 t2c = dvec2(trap2Vector.xy);

#define trap1Vector par[5]
  dvec2 trap1vector = dvec2(trap1Vector.xy);

  for( int i=0; i<iters; i++ ) {
    if( dot(z,z)>iters ) break;

    // Z' -> 2·Z·Z' + 1
    dz = 2.0*dvec2(z.x*dz.x-z.y*dz.y, z.x*dz.y + z.y*dz.x ) + dvec2(1.0,0.0);
      
    // Z -> Z² + c      
    z = cc + dvec2( z.x*z.x - z.y*z.y, 2.0*z.x*z.y );
      
    // trap 1
    double d1 = abs(dot(z-dvec2(0.0,1.0), trap1vector));
    double ff = step( d1, 1.0 );
    co2 += ff;
    trap1 += ff*d1;

    //trap2
    trap2 = min( trap2, dot(z-t2c,z-t2c) );
  }

    // distance, d(c) = |Z|·log|Z|/|Z'|
  double d = sqrt( dot(z,z)/dot(dz,dz) )*log(float(length(z)));

#define fzoo par[1].x  //{min=0 max=1 step=.001}
//  float zoo = 1.0 / 250.0;

  double zoo = dspeed * fzoo;
  
  float c1 = pow( float(clamp( 2.00*d/zoo,    0.0, 1.0 )), 0.5 );
  float c2 = pow( float(clamp( 1.5*trap1/co2, 0.0, 1.0 )), 2.0 );
  float c3 = pow( float(clamp( 0.4*trap2, 0.0, 1.0 )), 0.25 );

  vec3 col1 = 0.5 + 0.5*sin( 3.0 + 4.0*c2 + surfaceColor );
  vec3 col2 = 0.5 + 0.5*sin( 4.1 + 2.0*c3 + surfaceColor2 );
  vec3 col = 2.0*sqrt(c1*col1*col2);
  return col;
}
#endif


void main() {
  dvec3 p = deye;

  // Intersect the view ray with the 2D plane at Z==0
  double totalD = -p.z / dir.z;
  p += totalD * dir;

  vec3 col = vec3(0.0);
  if (totalD > 0.0) col = getColor2D(dvec2(p));

  // Write zBuffer and pixel
  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / float(clamp(totalD/length(dir), zNear, zFar)));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
