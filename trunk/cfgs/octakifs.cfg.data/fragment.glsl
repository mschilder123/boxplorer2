// Mandelbox shader by Rrrola
// Original mandelbox formula by Tglad
// Original kifs formulas by Knighty,parts by Marius
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal
// octakifs /bermarte/visual>goto fractalforums.com 
//#pragma OPENCL_EXTENSION cl_khr_fp64 : enable

#define P0 p0
#define Z0 z0
#define CX par[7].y// {min=-3.0 max=3 step=.001}
#define CY par[7].x// {min=-3.0 max=3 step=.001}
#define CZ par[8].x// {min=-3.0 max=3 step=.001}
#define Angle par[3].x// {min=-6.0 max=6 step=.001}
#define scale par[4].x// {min=-3 max=6 step=.001}
#define surfaceColor1 par[1]
#define surfaceColor2 par[2]
#define surfaceColor3 par[6]
#define glowColor par[5]
//later?
//#define fogStrength par[9].y  // {min=0 max=100 step=.25}
//vec3 z0 = normalize(vec3(par[8].y, par[3].y, par[3].x));
//#define LoDbase par[9].x  // {min=1 max=20 step=.25}
//#define LoDpow par[9].y  // {min=0 max=50 step=.25}
#define DE de_octa
#define COLOR color_octa
#define DIST_MULTIPLIER 1.0
#define MAX_DIST 4.0



// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[20];//4 later

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters,    // Number of fractal iterations.
  color_iters,        // Number of fractal iterations for coloring.
  max_steps,          // Maximum raymarching steps.
  frameno;

// Colors. Can be negative or >1 for interestiong effects.
//vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  //surfaceColor1 = vec3(0.95, 0.64, 0.1),
  //surfaceColor2 = vec3(0.89, 0.95, 0.75),
  //surfaceColor3 = vec3(0.55, 0.06, 0.03),
  //specularColor = vec3(1.0, 0.8, 0.4),
  //glowColor = vec3(0.03, 0.4, 0.4),
  //aoColor = vec3(.1, .1, .1);
//////////////
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  //surfaceColor1 = vec3(0.95, 0.12, 0.1),
  //surfaceColor2 = vec3(0.29, 0.4, 0.75),
  //surfaceColor3 = vec3(0.4, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  //glowColor = vec3(.1, 0.4, 0.4),
  aoColor = vec3(.1, .1, .1);

// precomputed constants
//float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
//vec4 scale = vec4(1.1, 1.1,1.1, abs(1.1)) / minRad2;
float absScalem1 = abs(1.1 - 1.0);//SCALE
float AbsScaleRaisedTo1mIters = pow(abs(1.1), float(1-45));
// Rotate around vector
//WAS float Angle = par[3].x;
float csat = cos(Angle);
float ssat = sin(Angle);
float usat = 1.0-cos(Angle);
vec3 z0 = normalize(vec3(par[8].y, par[3].y, par[3].x));
//WAS vec3 z0 = normalize(vec3(par[2].x, par[2].y, par[3].y));
mat3 RotationMatrix = mat3( z0.x*z0.x*usat + csat,      z0.x*z0.y*usat + z0.z*ssat, z0.x*z0.z*usat - z0.y*ssat,
                            z0.y*z0.x*usat - z0.z*ssat, z0.y*z0.y*usat + csat,      z0.y*z0.z*usat + z0.x*ssat,
			    z0.z*z0.x*usat + z0.y*ssat, z0.z*z0.y*usat - z0.x*ssat, z0.z*z0.z*usat + csat
		      );

// Compute the distance from `pos` to the Mandelbox.
float de_octa(vec3 z0) {
    float r=z0.x*z0.x+z0.y*z0.y+z0.z*z0.z;
    
    int i = 0;
    for (i=0;i< iters && r<20.0;i++){
      z0 *= RotationMatrix;
      vec3 zz0;
      if( z0.x+ z0.y<0.0){zz0.x=-z0.y;z0.y=-z0.x;z0.x=zz0.x;}
      if( z0.x+ z0.z<0.0){zz0.x=-z0.z;z0.z=-z0.x;z0.x=zz0.x;}
      if( z0.x- z0.y<0.0){zz0.x=z0.y;z0.y=z0.x;z0.x=zz0.x;}
      if( z0.x- z0.z<0.0){zz0.x=z0.z;z0.z=z0.x;z0.x=zz0.x;}
      //WAS CX=CY=CZ=1.
      z0.x=z0.x*scale-CX*(scale-CX);
      z0.y=z0.y*scale-CY*(scale-CY);
      z0.z=z0.z*scale-CZ*(scale-CZ);
      r=z0.x*z0.x+z0.y*z0.y+z0.z*z0.z;
      //z0 *= RotationMatrix;
}
  float k= (sqrt(r)-2.0)*pow(scale,float(-i));
  return max(k,-k);
}
// Compute the color at `pos`.
vec3 color_octa(vec3 pos) {
  vec3 p = pos, p0 = p;
  float trap = 1.0;
  for (int i=0; i<color_iters; i++) {
    p *= RotationMatrix;
    float r2 = dot(p, p);
    vec3 pp0;
    if( p.x+ p.y<0.0){pp0.x=-p.y;p.y=-p.x;p.x=pp0.x;}
    if( p.x+ p.z<0.0){pp0.x=-p.z;p.z=-p.x;p.x=pp0.x;}
    if( p.x- p.y<0.0){pp0.x=p.y;p.y=p.x;p.x=pp0.x;}
    if( p.x- p.z<0.0){pp0.x=p.z;p.z=p.x;p.x=pp0.x;}
    //WAS CX=CY=CZ=1.
    p.x=p.x*scale-CX*(scale-CX);
    p.y=p.y*scale-CY*(scale-CY);
    p.z=p.z*scale-CZ*(scale-CZ);
    r2=p.x*p.x+p.y*p.y+p.z*p.z;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
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
    -DE(pos-Eps.yxx)+DE(pos+Eps.yxx),
    -DE(pos-Eps.xyx)+DE(pos+Eps.xyx),
    -DE(pos-Eps.xxy)+DE(pos+Eps.xxy)

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
float ambient_occlusion(vec3 p, vec3 n) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = DE(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_eps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}


void main() {
  vec3 p = eye, dp = normalize(dir);

  float totalD = 0.0, D = 3.4e38, extraD = 0.0, lastD;

  // Intersect the view ray with the Mandelbox using raymarching.
  int steps;
  for (steps=0; steps<max_steps; steps++) {
    lastD = D;
    D = DE(p + totalD * dp);

    // Overstepping: have we jumped too far? Cancel last step.
    if (extraD > 0.0 && D < extraD) {
      totalD -= extraD;
      extraD = 0.0;
      D = 3.4e38;
      steps--;
      continue;
    }

    if (D < min_dist || D > MAX_DIST) break;

    totalD += D;

    // Overstepping is based on the optimal length of the last step.
    totalD += extraD = 0.096 * D*(D+extraD)/lastD;
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, D);
    col = COLOR(p);
    col = blinn_phong(n, -dp, normalize(eye+vec3(0,1,0)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > min_dist) {
      col = mix(col, backgroundColor, clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
