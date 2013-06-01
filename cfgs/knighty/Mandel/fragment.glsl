// 2.5D mandelbrot set (the lazy way ;o)) shader by knighty. Based on mandelbox shader by Rrrola.
// You could modify it to display julia sets 

/* Controls: numbers are the keys. x <=> left/right arrows. y <=> up/down arrows.
		Slope:       		par[0].x 
		Distance multiplier:    par[0].y 
	Have fun!
*/

#define DIST_MULTIPLIER par[0].y // {min=.001 max=1 step=.001}
#define DE_EPS 0.0001 
#define MAX_DIST 10.0

// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;           // {min=1e-6 max=.01 step=1e-6} Distance at which raymarching stops.
uniform float ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform float speed;

uniform int iters;    // {min=1 max=1000 step=1} Number of fractal iterations.
uniform int color_iters;        // Number of fractal iterations for coloring.
uniform int max_steps;          // Maximum raymarching steps.

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.03, 0.4, 0.4),
  aoColor = vec3(0, 0, 0);

//Fuctions to call.
#define d Mandel //d_PKlein//d_PZshape //
#define color color_0//color_PKlein//

#define BAILOUT 10000000.
#define SLOPE par[0].x // {min=-5 max=5 step=.01}
float Mandel(vec3 pos) {
	vec2 z=pos.xy;
	vec2 z0=z;
	float r2=dot(z,z);
#if 0
	vec2 dz=vec2(1.,0.);
	int i=0;
	for(i=0;i<iters && r2<BAILOUT;i++){
		dz=2.*vec2(dz.x*z.x-dz.y*z.y+1.,dz.x*z.y+dz.y*z.x);
		z=vec2(z.x*z.x-z.y*z.y,z.x*z.y*2.)+z0;
		r2=dot(z,z);
	}
	float r=sqrt(r2);
	float dr=length(dz);
	dr=DIST_MULTIPLIER*(r-2.)*log(r+1.)/dr;
	//if(r2<4.) dr=0.;
	dr=abs(dr); //max(dr,0.);
	return (pos.z-SLOPE*dr)/sqrt(1.+SLOPE*SLOPE);//transforms the 2D estimated dist to 3D
	//cone tracing would be faster but I'm too lazy
#else
	float dr=1.,ddr=0.,r;
	int i=0;
	for(i=0;i<iters && r2<BAILOUT;++i) {
		r=sqrt(r2);
		ddr = 2.*(dr*dr+r*ddr);
		dr = 2.*dr*r+1.;
		z=vec2(z.x*z.x-z.y*z.y,z.x*z.y*2.)+z0;
                r2 = dot(z,z);
	}
        r = sqrt(r2);
	float de = .5*log(dr)*(dr)/ddr;
	//float de = dr/ddr;
	return (pos.z-SLOPE*de)/sqrt(1.+SLOPE*SLOPE);
#endif
}
// Compute the color at `pos`. I'm too lazy (once again lol) to do a good one!
vec3 color_0(vec3 pos) {
	return vec3(0.7,0.6,0.4);
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


// FAKE Ambient occlusion approximation.
// uses current distance estimate as first dist. the size of AO is independent from distance from eye
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side) {
  float ao_ed=DistAtp*ao_eps/min_dist;//Dividing by min_dist makes the AO effect independent from changing min_dist
  float ao = 1.0, w = ao_strength/ao_ed;//ps;
  float dist = 2.0 * ao_ed;//ps;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 1.0/2.0;
    dist = dist*2.0 - ao_ed;//ps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

#define ONE_PLUS_ULP 1.000000059604644775390625 //1.00000011920928955078125 // 
#define ONE_MINUS_ULP 0.999999940395355224609375 //0.99999988079071044921875 // 
#define ULP 0.000000059604644775390625 //0.00000011920928955078125 // 
void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);
  vec3 p = eye_in, dp = normalize(dir);

  float totalD = 0.0, D = 1000.0, extraD = 0.0, lastD;
  D = d(p + totalD * dp);
  float side = sign(D);
  D = abs(D); 
  float MINDIST_MULT=1.0/(1.0+min_dist);
  totalD = D *0.5* MINDIST_MULT;

  // Intersect the view ray with the Mandelbox using raymarching.
  // The distance field actually marched is the "calculated DE" minus (totalD * min_dist)
  // A perfect distance field have a gradient magnitude = 1. Assuming d() gives a perfect DE, 
  //we have to multiply D with MINDIST_MULT in order to restore a gradient magnitude of 1
  int steps;

  for (steps=0; steps<max_steps && abs(D) > min_dist*totalD*1.0/8.0 && abs(D)>256.0*ULP && totalD < MAX_DIST; steps++) {
    D = (side * d(p + totalD * dp) - totalD * min_dist) * MINDIST_MULT;
    totalD+=D;
  }
  p += totalD * dp;
  //normal_eps = max(normal_eps,1.0*D);
  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    //D = max(D, max(DE_EPS*totalD,0.0000152587890625));
	float D1 = min_dist*1.0/32.0*totalD;
	vec3 n = side * normal(p, D1);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(/*eye_in+*/vec3(1.0,0.5,0.7)/*+dp*/), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, D1, side));

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
