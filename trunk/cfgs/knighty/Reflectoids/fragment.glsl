// "Pseudo Kleinian" shader by knighty. Based on mandelbox shader by Rrrola.
// Original formula by Theli-at. Actally "reverse ingeneered" from his MB3D param files. Thanks theli :)
// - http://www.fractalforums.com/mandelbulb-3d/mandelbulb-3d-parameters-list/msg28870/#msg28870

/* Controls: numbers are the keys. x <=> left/right arrows. y <=> up/down arrows.
		Size of the box folding cell:   vec3(par[1].y,par[1].x,par[2].y)
		Size of inversion sphere:       par[0].y //Well it is redundant with the previous but I think it's convenient
		Julia seed:                     vec3(par[2].x,par[3].y,par[3].x)
		Translation of the basic shape: vec3(par[4].y,par[4].x,par[5].y)
		Chage a little the basic shape: par[5].x
		A small value added to the DE:  par[0].x
		Color parameter:                vec3(par[6].y,par[6].x,par[7].y)
		Blend with 0 trap:              par[7].x
	Have fun!
*/

#define DIST_MULTIPLIER par[9].x  // {min=.01 max=1 step=.01}
#define MAX_DIST 10.0

#define INOUT(a,b) inout a b

#include "setup.inc"
#line 24

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;           // Distance at which raymarching stops.
uniform float ao_eps;             // Base distance at which ambient occlusion is estimated.
uniform float ao_strength;        // Strength of ambient occlusion.
uniform float glow_strength;      // How much glow is applied after max_steps.
uniform float dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters;                // Number of fractal iterations.
uniform int color_iters;          // Number of fractal iterations for coloring.
uniform int max_steps;            // Maximum raymarching steps.

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(0.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.03, 0.4, 0.4),
  aoColor = vec3(0, 0, 0);

//Fuctions to call.
#define d d_PKlein//d_PZshape //
#define d2 d_PKlein2//d_PZshape //
#define color color_PKlein//color_0//

// Compute the distance from `pos` to the PKlein basic shape.
float d_PZshape(vec3 p) {
#define ZHole par[9].y  //{min=-5 max=5 step=.01}
   float rxy = sign(ZHole)*(length(p.xy)-abs(ZHole));

#define TThickness par[5].x //{min=-5 max=5 step=.01}
#define Zmult par[8].y //{min=-5 max=5 step=.01}
#define Ziter par[8].x  //{min=1 max=100 step=1}

   for(int i=0; i<int(Ziter); i++) p.z=2.*clamp(p.z, -Zmult, Zmult)-p.z;

   // abs() causes banding at times.
   return max(rxy, abs(length(p.xy)*p.z-TThickness) / sqrt(dot(p,p)+abs(TThickness)));
}

// Compute the distance from `pos` to the PKlein.
float d_PKlein(vec3 p) {
   //Just scale=1 Julia box
	float r2=dot(p,p);
	float DEfactor=1.;

#define CSize vec3(par[1].y,par[1].x,par[2].y)
#define Size par[0].y //{min=-5 max=5 step=.01}
#define C vec3(par[2].x,par[3].y,par[3].x)
#define Offset vec3(par[4].y,par[4].x,par[5].y)
#define DEoffset par[0].x //{min=-5 max=5 step=.01}

	for(int i=0;i<iters;i++){
		//Box folding
		p=2.*clamp(p, -CSize, CSize)-p;
		//Inversion
		r2=dot(p,p);
		float k=max(Size/r2,1.);
		p*=k;DEfactor*=k;
		//julia seed
		p+=C;
	}
	//Call basic shape and scale its DE
	//Ok... not perfect because the inversion transformation is tricky
	//but works Ok with this shape (maybe because of the "tube" along Z-axis
	//You may need to adjust DIST_MULTIPLIER (par[9].x) especialy with non zero Julia seed
	return
	  (DIST_MULTIPLIER*d_PZshape(p-Offset)/abs(DEfactor)-DEoffset);
}

// Compute the distance from `pos` to the PKlein basic shape.
float d_PZshape2(vec3 p) {
   float rxy = sign(ZHole)*(length(p.xy)-abs(ZHole));

   int i = 0;
   for(; i<int(Ziter); i++) p.z=2.*clamp(p.z, -Zmult, Zmult)-p.z;
   //i = 2*(i&1) - 1;

   // abs() causes banding at times.
   return max(rxy, abs/*-float(i)* */(length(p.xy)*p.z-TThickness) / sqrt(dot(p,p)+abs(TThickness)));
}

// Compute the distance from `pos` to the PKlein.
float d_PKlein2(vec3 p) {
   //Just scale=1 Julia box
	float r2=dot(p,p);
	float DEfactor=1.;

	for(int i=0;i<iters;i++){
		//Box folding
		p=2.*clamp(p, -CSize, CSize)-p;
		//Inversion
		r2=dot(p,p);
		float k=max(Size/r2,1.);
		p*=k;DEfactor*=k;
		//julia seed
		p+=C;
	}
	//Call basic shape and scale its DE
	//Ok... not perfect because the inversion transformation is tricky
	//but works Ok with this shape (maybe because of the "tube" along Z-axis
	//You may need to adjust DIST_MULTIPLIER (par[9].x) especialy with non zero Julia seed
	return
	  (DIST_MULTIPLIER*d_PZshape2(p-Offset)/abs(DEfactor)-DEoffset);
}

// Compute the color.
vec3 color_PKlein(vec3 p) {
   //Just scale=1 Julia box
	float r2=dot(p,p);
	float DEfactor=1.;
	vec4  col=vec4(0.0);
	float rmin=10000.0;//r2;
/*#define CSize vec3(par[1].y,par[1].x,par[2].y)
#define Size par[0].y
#define C vec3(par[2].x,par[3].y,par[3].x)
#define Offset vec3(par[4].y,par[4].x,par[5].y)
#define DEoffset par[0].x*/
	for(int i=0;i<color_iters /*&& r2<1000000.0*/;i++){
		//Box folding
		vec3 p1=2.*clamp(p, -CSize, CSize)-p;
		col.xyz+=abs(p-p1);//vec3(notEqual(p,p1));
		p=p1;
		//Inversion
		r2=dot(p,p);
		float k=max(Size/r2,1.);
		col.w+=abs(k-1.);
		p*=k;DEfactor*=k;
		//julia seed
		p+=C;
		r2=dot(p,p);
		rmin=min(rmin,r2);
	}
	//rmin=min(1.,r2);
	return mix(vec3(sqrt(rmin)),(0.5+0.5*sin(col.z*vec3(par[6].y,par[6].x,par[7].y))),par[7].x);//vec3(sqrt(rmin));//*col.xyz/(iters+1.);
}
// Compute the color at `pos`. I'm too lazy to do a good one!
vec3 color_0(vec3 pos) {
	return vec3(0.7,0.6,0.4);
}

float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  //vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  vec2 Eps = vec2(0, max(d_pos, normal_eps));
  return normalize(vec3(
    -d2(pos-Eps.yxx)+d2(pos+Eps.yxx),
    -d2(pos-Eps.xyx)+d2(pos+Eps.xyx),
    -d2(pos-Eps.xxy)+d2(pos+Eps.xxy)
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
float marche(inout vec3 p, in vec3 dp, inout float D, inout float totalD, in float side, in float MINDIST_MULT){
	// Intersect the view ray with the Mandelbox using raymarching.
	// The distance field actually marched is the "calculated DE" minus (totalD * min_dist)
	// A perfect distance field have a gradient magnitude = 1. Assuming d() gives a perfect DE, 
	//we have to multiply D with MINDIST_MULT in order to restore a gradient magnitude of 1
	int steps;
	for (steps=0; steps<max_steps && abs(D)>max(totalD*8192.0*ULP,ULP) && totalD < MAX_DIST; steps++) {
		totalD+=D;
		D = (side * d(p + totalD * dp) - totalD * min_dist) * MINDIST_MULT;
	}
	p += (totalD+D) * dp;
	return float(steps);
}

float hash( float n ) {
    return fract(sin(n)*5345.8621276);
}

float noise( in vec2 x ) {
    vec2 p = floor(x);
    vec2 f = fract(x);

    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*61.0;

    float res = mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                    mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y);

    return fract(res);
}

void main() {
  vec3 eye_in, dp; 

  if (!setup_ray(eye, dir, eye_in, dp)) return;

  float noise = noise(gl_FragCoord.xy / vec2(xres, yres));

  vec3 p = eye_in; // + dp * min_dist * 0.75;

  float totalD = 0.0, D = d(p);
  float side = sign(D);
  D = noise * abs(D); 
  float MINDIST_MULT=1.0/(1.0+min_dist);
  D *= MINDIST_MULT;

	vec3 finalcol=vec3(0.);
	float refpart=1.0;
#define REFACTOR par[8].z  //{min=0 max=1 step=.01}
#define REFITER par[9].z
	int i=0;
	bool cont=true;
     float firstD = 0.;  // save first step for depth buffer
while(i<int(REFITER) && cont){
  float steps=marche(p, dp, D, totalD, side, MINDIST_MULT);
  if (i == 0) firstD = totalD + D;
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (totalD < MAX_DIST) {
	float D1 = min_dist*.5*totalD;
	vec3 n = side * normal(p, max(256.0*ULP,D1));
    col = color(p);
    col = blinn_phong(n, -dp, normalize(/*eye+*/vec3(1.0,0.5,0.7)/*+dp*/), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, D1, side));
	
	// update the ray
	dp=reflect(dp,n);

	p -= (totalD+D) * dp;
//	D=2.*D1;
	D = (9. + noise) * D1;

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > max(totalD*8192.0*ULP,ULP)){//min_dist) {
      col = mix(col, backgroundColor, clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
    }
	//col = mix(col, glowColor, steps/float(max_steps) * glow_strength);
  }else {
	cont=false;
	}

  // Glow is based on the number of steps.
  col = mix(col, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);
  finalcol+=refpart*col;
  refpart*=REFACTOR;
  i++;
}

  write_pixel(dir, firstD, finalcol);
}
