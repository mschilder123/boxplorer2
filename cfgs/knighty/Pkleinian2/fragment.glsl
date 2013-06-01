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
#define d d_PKlein//d_PZshape //
#define color color_PKlein//color_0//

// Compute the distance from `pos` to the PKlein basic shape.
float d_PZshape(vec3 p) {
   float rxy=sign(par[9].y)*(length(p.xy)-abs(par[9].y));
#define TThickness par[5].x
#define Zmult par[8].y
#define Ziter int(par[8].x)
   for(int i=0; i<Ziter; i++) p.z=2.*clamp(p.z, -Zmult, Zmult)-p.z;
   return max(rxy,abs(length(p.xy)*p.z-TThickness) / sqrt(dot(p,p)+abs(TThickness)));
}

// Compute the distance from `pos` to the PKlein.
float d_PKlein(vec3 p) {
   //Just scale=1 Julia box
	float r2=dot(p,p);
	float DEfactor=1.;
#define CSize vec3(par[1].y,par[1].x,par[2].y)
#define Size par[0].y
#define C vec3(par[2].x,par[3].y,par[3].x)
#define Offset vec3(par[4].y,par[4].x,par[5].y)
#define DEoffset par[0].x
	for(int i=0;i<iters /*&& r2<1000000.0*/;i++){
		//Box folding
		p=2.*clamp(p, -CSize, CSize)-p;
		//Inversion
		r2=dot(p,p);
		float k=max(Size/r2,1.);
		p*=k;DEfactor*=k;
		//julia seed
		p+=C;
		//r2=dot(p,p);
	}
	//Call basic shape and scale its DE
	//Ok... not perfect because the inversion transformation is tricky
	//but works Ok with this shape (maybe because of the "tube" along Z-axis
	//You may need to adjust DIST_MULTIPLIER (par[9].x) especialy with non zero Julia seed
	return (DIST_MULTIPLIER*d_PZshape(p-Offset)/DEfactor-DEoffset);
}

// Compute the color.
vec3 color_PKlein(vec3 p) {
   //Just scale=1 Julia box
	float r2=dot(p,p);
	float DEfactor=1.;
	vec4  col=vec4(0.0);
	float rmin=r2;
/*#define CSize vec3(par[1].y,par[1].x,par[2].y)
#define Size par[0].y
#define C vec3(par[2].x,par[3].y,par[3].x)
#define Offset vec3(par[4].y,par[4].x,par[5].y)
#define DEoffset par[0].x*/
	for(int i=0;i<iters /*&& r2<1000000.0*/;i++){
		//Box folding
		vec3 p1=2.*clamp(p, -CSize, CSize)-p;
		col.xyz+=abs(p-p1);//vec3(notEqual(p,p1));
		p=p1;
		//Inversion
		r2=dot(p,p);
		float k=max(Size/r2,1.);col.w+=abs(k-1.);
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
