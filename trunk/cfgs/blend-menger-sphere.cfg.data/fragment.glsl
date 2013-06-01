// Mandelbox shader by Rrrola. Modified by knighty :o)
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal

#define P0 p0                    // standard Mandelbox
//#define P0 vec4(par[1].x,par[1].y,par[2].y,1)  // Mandelbox Julia

#define SCALE par[0].y
#define MINRAD2 par[0].x

#define MENGER_SIZ par[1].y  //{min=0 max=10 step=.01}
#define MENGER_SX par[2].y  //{min=0 max=10 step=.01}
#define MENGER_SY par[2].x  //{min=0 max=10 step=.01}
#define MENGER_SZ par[3].y  //{min=0 max=10 step=.01}

#define DIST_MULTIPLIER par[9].x //{min=0 max=10 step=.01}
#define DE_EPS 0.0001 //par[9].y
#define MAX_DIST 10.0

// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform float speed;

uniform int iters,    // Number of fractal iterations.
  color_iters,        // Number of fractal iterations for coloring.
  max_steps;          // Maximum raymarching steps.

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.03, 0.4, 0.4),
  aoColor = vec3(0, 0, 0);

// precomputed constants
float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
vec4 scale = vec4(SCALE, SCALE, SCALE, abs(SCALE)) / minRad2;
float absScalem1 = abs(SCALE - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(SCALE), float(1-iters));

//Fuctions to call.
#define d d_MengerInterMandelBox //d_Mengersphere //d_Menger2 //d_mandelBox//d_MengerInterMandelBox//d_MengerIntersphere //
#define color color_mandelBox//

// Compute the distance from `pos` to the sphere.
float d_sphere(vec3 pos) {
   return length(pos)-par[4].x;
}
// Compute the distance from `pos` to the Menger.
// we carve at each iteration a 3D "+" shapes which scale depends on the iteration number
// The mod() allow us to make multiple copies of the caved shape
float d_Menger2(vec3 pos) {//Old one
   pos=pos*0.5+vec3(0.5);
   vec3 p=abs(pos-vec3(0.5))-vec3(0.5);
   float di1=max(p.x,max(p.y,p.z));
   float di=di1;//0.0;//
   float pf=1.0,pfi=1.0,dmax=1.0/6.0;
   for(int m=1;m<iters && dmax>di;m++){
    vec3 pa=mod(3.0*pf*pos,3.0);
      pf*=3.0;pfi*=1.0/3.0;dmax*=1.0/3.0;
      p=0.5-abs(pa+vec3(-1.5));
      p=max(p.xyz,p.yzx);
      di1=min(p.x,min(p.y,p.z))*pfi;
      di=max(di,di1);
   }
   return di*2.0;
}
//Here I've added translation (2 and 3 keys) to pos and "+" shape size modifier (0 key).
// One can also add rotations.
// also there are much oter possibilities:
// - other shapes than the 3D "+" (see buddhi's spherical menger in "Revenge of half eaten menger" thread at www.fractalforums.com
// - size, translation, rotation and carved shape may be different at each iteration.
// - if possible, other "foldings" than mod() (which gives cubic repetition)
float d_Menger(vec3 pos) {//modified one
   pos=pos*0.5+vec3(0.5);
   vec3 p=abs(pos-vec3(0.5))-vec3(0.5);
   float di1=2.0*max(p.x,max(p.y,p.z));
   float di=di1;//0.0;//
   float pf=3.0,pfi=1.0/3.0,k=MENGER_SIZ;
   for(int m=1;m<iters && di<k*pfi;m++){
    pos += vec3(MENGER_SX,MENGER_SY,MENGER_SZ); 
      vec3 pa=mod(pf*pos,3.0);
      p=k-abs(2.0*pa+vec3(-3.0));
      p=max(p.xyz,p.yzx);
      di1=min(p.x,min(p.y,p.z))*pfi;
      di=max(di,di1);
    pf*=3.0;pfi*=1.0/3.0;
   }
   return di;
}
// Compute the color at `pos`.
vec3 color_Menger(vec3 pos) {
  return vec3(0.7,0.6,0.4);
}

// Compute the distance from `pos` to the mix of Menger sponge and sphere.
float d_Mengersphere(vec3 pos) {
   return mix(d_Menger(pos),d_sphere(pos),par[4].y)/max(1.0,2.0*abs(par[4].y-0.5));
}


// Compute the distance from `pos` to the Mandelbulb.
float d_mandelBulb(vec3 pos) {
  return 0.0;
}
// Compute the color at `pos`.
vec3 color_mandelBulb(vec3 pos) {
  return vec3(0.0,0.0,0.0);
}
// Compute the distance from `pos` to the Mandelbox.
float d_mandelBox(vec3 pos) {
  vec4 p = vec4(pos,1), p0 = p;  // p.w is the distance estimate
  float r2 = dot(p.xyz, p.xyz);
  int i=0; 
  for (i=0; i<iters && r2<10000.0; i++) {
    // box folding: if (p>1) p = 2-p; else if (p<-1) p = -2-p;
//    p.xyz = abs(1.0+p.xyz) - p.xyz - abs(1.0-p.xyz);  // add;add;abs.add;abs.add (130.4%)
//    p.xyz = clamp(p.xyz*0.5+0.5, 0.0, 1.0) * 4.0 - 2.0 - p.xyz;  // mad.sat;mad;add (102.3%)
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;  // min;max;mad

    // sphere folding: if (r2 < minRad2) p /= minRad2; else if (r2 < 1.0) p /= r2;
    r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);  // dp3,div,max.sat,mul

    // scale, translate
    p = p*scale + P0;
    r2 = dot(p.xyz, p.xyz);
  }
  return ((length(p.xyz) - absScalem1) / p.w /*- pow(abs(SCALE), float(1-i))/*AbsScaleRaisedTo1mIters*/) * DIST_MULTIPLIER;
}


// Compute the color at `pos`.
vec3 color_mandelBox(vec3 pos) {
  vec3 p = pos, p0 = p;
  float trap = 100.0;

  for (int i=0; i<color_iters; i++) {
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale.xyz + P0.xyz;
    r2 = dot(p.xyz, p.xyz);
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, 0.25*sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}

// Compute the distance from `pos` to the intersection of Menger sponge and Mandelbox.
float d_MengerInterMandelBox(vec3 pos) {
   return min(d_Menger(pos),d_mandelBox(pos));
}

// Compute the distance from `pos` to the intersection of Menger sponge and sphere.
float d_MengerIntersphere(vec3 pos) {
   return min(d_mandelBox(pos),d_sphere(pos));
}


float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  return normalize(vec3(
  // 2-tap forward differences, error = O(eps)
//    -d_pos+d(pos+Eps.yxx),
//    -d_pos+d(pos+Eps.xyx),
//    -d_pos+d(pos+Eps.xxy)

  // 3-tap central differences, error = O(eps^2)
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)

  // 4-tap forward differences, error = O(eps^3)
//    -2.0*d(pos-Eps.yxx)-3.0*d_pos+6.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
//    -2.0*d(pos-Eps.xyx)-3.0*d_pos+6.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
//    -2.0*d(pos-Eps.xxy)-3.0*d_pos+6.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)

  // 5-tap central differences, error = O(eps^4)
//    d(pos-Eps.zxx)-8.0*d(pos-Eps.yxx)+8.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
//    d(pos-Eps.xzx)-8.0*d(pos-Eps.xyx)+8.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
//    d(pos-Eps.xxz)-8.0*d(pos-Eps.xxy)+8.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)
  ));
}


// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
}


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side) {
  float ao_ed=DistAtp*100.0;
  float ao = 1.0, w = ao_strength/ao_ed;//ps;
  float dist = 2.0 * ao_ed;//ps;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;//ps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 30x speed.
void setup_stereo(inout vec3 eye_in, inout vec3 dp) {
#if !defined(ST_NONE)
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 4.0 * (fract(gl_FragCoord.y * 0.5) - .5) * abs(speed), 0, 0, 0));
#else
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0, 0, 0));
#endif
  eye_in = eye + eye_d;
  dp = normalize(dir * (focus + 30.0) * abs(speed) - eye_d);
#else  // ST_NONE
  eye_in = eye;
  dp = normalize(dir);
#endif
}

void main() {
  vec3 eye_in, dp; setup_stereo(eye_in, dp);
  vec3 p = eye_in;

  float totalD = 0.0, D = 1000.0, extraD = 0.0, lastD;
  D = d(p + totalD * dp);
  float side = sign(D);
  D = D*side; D *= 0.95;
  totalD = D;

  // Intersect the view ray with the Mandelbox using raymarching.
  int steps;
  for (steps=0; steps<max_steps && D > max(DE_EPS*totalD,0.0000152587890625) && totalD < MAX_DIST; steps++) {
    D = side * d(p + totalD * dp);
    D = D - min_dist * (totalD + D) ;
    totalD += D;
  }

  p += totalD * dp;
  //normal_eps = max(normal_eps,1.0*D);
  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    D = max(D, max(DE_EPS*totalD,0.0000152587890625));
    vec3 n = side * normal(p, D);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(/*eye_in+*/vec3(1.0,0.5,0.7)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, D, side));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > min_dist) {
      col = mix(col, backgroundColor, clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  float zNear = abs(speed);
  float zFar = 65535.0 * zNear;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
