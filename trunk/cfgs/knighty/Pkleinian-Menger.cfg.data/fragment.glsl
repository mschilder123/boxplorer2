// "Pseudo Kleinian" with Menger K-IFS shader by knighty. Based on mandelbox shader by Rrrola.
// Original formula by Theli-at. Actally "reverse ingeneered" from his MB3D param files. Thanks theli :)
// - http://www.fractalforums.com/mandelbulb-3d/mandelbulb-3d-parameters-list/msg28870/#msg28870

/* Controls: numbers are the keys. x <=> left/right arrows. y <=> up/down arrows.
	  Pseudo Kleinian:	
		Size of the box folding cell:   vec3(par[1].y,par[1].x,par[2].y)
		Size of inversion sphere:       par[0].y //Well it is redundant with the previous but I think it's convenient
		Julia seed:                     vec3(par[2].x,par[3].y,par[3].x)
		Translation of the basic shape: vec3(par[4].y,par[4].x,par[5].y)
		A small value added to the DE:  par[0].x
	  KIFS Mengers:
		[5].x Size of the Menger
		[6].y x component of center of scaling point
		[6].x y component of center of scaling point
		[7].y z component of center of scaling point
		[7].x 
		[8].y 
		[8].x 
		[9].y Iterations number
		[9].x Distance Multiplier 
		[0].x scale
	Have fun!
*/

#define DIST_MULTIPLIER par[9].x // {min=.001 max=1 step=.001}
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

//Fuctions to call.
#define d d_PKlein//de_KIFSMenger //
#define color color_0

//KIFS Mengers:
//controls:
//[6].y x component of center of scaling point
//[6].x y component of center of scaling point
//[7].y z component of center of scaling point
//[7].x x component of one more box folding 
//[8].y y component of one more box folding 
//[8].x z component of one more box folding 
//[9].y Iterations number
//[0].x scale
float de_KIFSMenger(vec3 z0) {// kaleidoscopic IFS Menger box
#define MenIters par[9].y
#define MenSCALE par[0].x
#define Siz par[5].x//0.25//
//#define CSiz vec3(par[7].x,par[8].y,par[8].x)
	//z0=2.*clamp(z0, -CSiz, CSiz)-z0;//Got compilation exception :/, maybe out of registers?
	z0*=1.0/Siz;
	int i;
	vec3 kz=abs(z0);
	float r=max(kz.x,max(kz.y,kz.z));
	for (i=0;i<int(MenIters) && r<100.0;i++){
		z0=abs(z0);
		if( z0.x - z0.y < 0.0) z0.xy=z0.yx;
		if( z0.x - z0.z < 0.0) z0.xz=z0.zx;//for full symmetry
		if( z0.y - z0.z < 0.0) z0.yz=z0.zy;
		if( z0.z - 1.0/3.0 < 0.0) z0.z += 2.0/3.0 - 2.0 * z0.z;
		z0 = z0 * MenSCALE - vec3(par[6].y,par[6].x,par[7].y) * vec3(MenSCALE-1.0);
		kz=abs(z0);
		r=max(kz.x,max(kz.y,kz.z));
	}
	kz=abs(vec3(par[6].y,par[6].x,par[7].y));
	float sub=max(kz.x,max(kz.y,kz.z));
	return (r-sub)*pow(MenSCALE,-float(i))*Siz;
}

// Compute the distance from `pos` to the PKlein.
float d_PKlein(vec3 p) {
   //Just scale=1 Julia box
	float DEfactor=1.;
#define CSize vec3(par[1].y,par[1].x,par[2].y)
#define Size par[0].y
#define C vec3(par[2].x,par[3].y,par[3].x)
#define Offset vec3(par[4].y,par[4].x,par[5].y)
	for(int i=0;i<iters;i++){
		//Box folding
		p=2.*clamp(p, -CSize, CSize)-p;
		//Sphere folding
		float r2=dot(p,p);
		float k=max(Size/r2,1.);
		p*=k;DEfactor*=k;
		//julia Seed
		p+=C;
	}
	//Call basic shape and scale its DE
	return abs(DIST_MULTIPLIER*de_KIFSMenger(p-Offset)/DEfactor);
}

// Compute the color at `pos`.
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
