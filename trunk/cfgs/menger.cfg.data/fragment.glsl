// menger and then some shader.
// Original shader by rrrola for mandelbox
// bermarte, marius: formulae from Knighty
// marius: refactored w/ reflections, background dome, ssponge, combi etc.
// marius: massaged so can be compiled as plain C++.

#ifndef _FAKE_GLSL_

#extension GL_ARB_shader_texture_lod : enable

#define DECLARE_DE(x)
#define DECLARE_COLORING(x)

#define INOUT(a,b) inout a b
#define IN(a,b) in a b
#define OUT(a,b) out a b

// distance estimator func
#ifndef d
#define d de_menger // PKlein,combi,menger,mandelbox,ssponge
#endif

#define DE_FUNC_VEC3 de_menger

// surface coloring func
#ifndef c
#define c c_menger  // PKlein,menger
#endif

#endif  // _FAKE_GLSL_

#define MAX_DIST 10.0
#define ULP 0.000000059604644775390625

#ifndef PI
#define PI 3.14159265
#endif

#include "setup.inc"
#line 41

// Camera position and direction.
//varying vec3 eye;
//varying vec3 dir;

//uniform float xres, yres, time, speed;

// Interactive parameters.
uniform vec3 par[20];

uniform float min_dist;           // Distance at which raymarching stops.
uniform float ao_eps;             // Base distance at which ambient occlusion is estimated.
uniform float ao_strength;        // Strength of ambient occlusion.
uniform float glow_strength;      // How much glow is applied after max_steps.
uniform float dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters;        // {min=1 max=1000} Number of fractal iterations.
uniform int color_iters;  // {min=0 max=1000} Coloration iterations.
uniform int max_steps;    // {min=1 max=1000} Maximum raymarching steps.
uniform int nrays;        // {min=1 max=10} # of ray bounces.

uniform bool julia;
#define JuliaVector par[1]

// Colors. Can be negative or >1 for interesting effects.
#define backgroundColor par[10]
#define surfaceColor1 par[11]
#define surfaceColor2 par[12]
#define surfaceColor3 par[13]
#define specularColor par[14]
#define glowColor par[15]
#define aoColor par[16]

#define DIST_MULTIPLIER par[8].z // {min=0.001 max=1 step=.001}

// Reflection params
#define SHINE par[9].x  // {min=0 max=1 step=.01} reflection contribution drop-off
#define SIGMA par[9].y  //{min=0 max=1 step=.01} surface roughness
#define ALBEDO par[9].z  //{min=0 max=1 step=.01} surface reflection

// Background blur
#define BG_BLUR par[1].z  // {min=0 max=8 step=.1}

// SphereSponge params
#define SS_MOD par[8].x  // {min=1 max=5 step=.01}
#define SS_SCALE par[8].y  // {min=1 max=5 step=.01}
#define SS_INTRA_D par[8].z  // {min=.01 max=10.0 step=.01}
#define SS_INITIAL_K par[0].z  // {min=0 max=5 step=.01}

// Mandelbox params
#define MB_SCALE par[0].y  //{min=-3 max=3 step=.01}
#define MB_MINRAD2 par[0].x  //{min=0 max=1 step=.001}

// Menger params
#define ME_offsetVector par[2]
#define ME_ITERS par[1].x // {min=0 max=20 step=1} iterations
#define ME_SIZE par[1].y // {min=0 max=100 step=.01} size of menger box

// PKleinian params
#define PZ_thickness par[5].x // {min=0.0 max=3 step=.001}
#define PZ_mult par[8].y // {min=0.0 max=1 step=.001}
#define PZ_iter par[8].x  // {min=0 max=10 step=1}
#define PZ_rxy par[7].y  // {min=0 max=25 step=.01}

#define PK_BSize vec3(par[1].y,par[1].x,par[2].y)
#define PK_CSize par[0].y // {min=0.0 max=10.0 step=.001}
#define PK_CVector par[3]
#define PK_Offset vec3(par[4].y,par[4].x,par[5].y)
#define PK_Color par[17]
#define PK_CMix par[7].x  // {min=-2 max=2 step=.01}
#define PK_DEoffset par[0].x

#define SPHERE_SZ par[4].z // {min=0 max=5 step=.01}

#define CYL_Vector par[17]

#define L1_Vector par[19]
#define L1_Size par[18].x  // {min=0.001 max=1.0 step=.001}

#if 0
float de_menger(vec3 z0) {
   float alpha = ME_SIZE;
   float x = z0.x;
   float y = z0.y;
   float z = z0.z;
   float scl=(alpha+2);
   float scl1=scl/alpha;
   float a=1/scl1;
   float r=x*x+y*y+z*z,dd=1.0;
   for(int i=0;i<int(ME_ITERS) && r<100;i++){
      x=abs(x);y=abs(y);z=abs(z);
      if(x<y){float t=x;x=y;y=t;}
      if(z<x){float t=z;z=x;x=t;}
      if(x<y){float t=x;x=y;y=t;}

      if(y<a && x>1-3*a+y){
         x-=1;z-=1;
         x*=scl1;y*=scl1;z*=scl1;dd*=scl1;
         x+=1;z+=1;
      }else{
         x-=1;y-=1;z-=1;
         x*=scl;y*=scl;z*=scl;dd*=scl;
         x+=1;y+=1;z+=1;
      }
      r=x*x+y*y+z*z;
   }
   return (sqrt(r)-1.75)/dd;
}
#else
float de_menger(vec3 z0) {
  float scale=2.;  // size of holes
  float r=1.0;
  z0 /= ME_SIZE;
  int n_iters = int(ME_ITERS);
  for (int i=0;i<n_iters;i++) {
    z0 = abs(z0) + ME_offsetVector;
    if( z0.x < z0.y){z0 = z0.yxz;}
    if( z0.x < z0.z){z0 = z0.zyx;}
    if( z0.y < z0.z){z0 = z0.xzy;}
    r = z0.x;
    z0 += z0*scale - scale;
    if(z0.z<-0.5*scale) z0.z+=scale;
  }
  return (r-1.0)*pow(float(scale+1.),float(1-n_iters))*ME_SIZE;
}
#endif
DECLARE_DE(de_menger)

//2D coordinates infinite tiling
vec2 Tile2D(vec2 p, vec2 a){
	p+=a;
	p-=a*floor(p/(a*4.0))*4.0;
	p-=max(p-a*2.0,0.0)*2.0;
	p-=a;
	return p;
}

float de_KIFSMenger(vec3 z0) {
	z0 /= ME_SIZE;
	z0=vec3(Tile2D(z0.xy,vec2(par[8].x, par[8].y)), z0.z); 
	int i;
	vec3 kz=abs(z0);
	const float ME_SCALE = 2.0;
	float r=max(kz.x,max(kz.y,kz.z));
	int n_iters = int(ME_ITERS);
	for (i=0;i<n_iters && r<10.0;i++){
		z0=abs(z0);
		if( z0.x - z0.y < 0.0) z0=z0.yxz;
		if( z0.x - z0.z < 0.0) z0=z0.zyx;
		if( z0.y - z0.z < 0.0) z0=z0.xzy;
		if( z0.z - 1.0/3.0 < 0.0) z0.z += 2.0/3.0 - 2.0 * z0.z;
		z0 += z0*ME_SCALE - ME_SCALE;
		kz=abs(z0);
		r=max(kz.x,max(kz.y,kz.z));
	}
	return (r-1.0)*pow(float(ME_SCALE+1.0),-float(i))*ME_SIZE;
}

vec3 c_menger(vec3 p) {
  return surfaceColor1;  // Boring but reflections make it interesting.
}
DECLARE_COLORING(c_menger)

// Compute the distance from `pos` to the PKlein basic shape.
float de_PZshape(vec3 p) {
  float rxy=sign(PZ_rxy)*(length(p.xy)-abs(PZ_rxy));
  for(int i=0; i<int(PZ_iter); i++) p.z=2.*clamp(p.z, -PZ_mult, PZ_mult)-p.z;
  return max(rxy,abs(length(p.xy)*p.z-PZ_thickness) / sqrt(dot(p,p)+abs(PZ_thickness)));
}

// Compute the distance from `pos` to the PKlein.
float de_PKlein(vec3 p) {
  //Just scale=1 Julia box
  float r2=dot(p,p);
  float DEfactor=1.;

  for(int i=0;i<iters;i++) {
    p=clamp(p, -PK_BSize, PK_BSize)*2.-p;  // Box folding
    r2=dot(p,p);  // Inversion
    float k=max(PK_CSize/r2,1.);
    p*=k; DEfactor*=k;
    p+=PK_CVector;  // julia seed
  }
  //Call basic shape and scale its DE
  //Ok... not perfect because the inversion transformation is tricky
  //but works Ok with this shape (maybe because of the "tube" along Z-axis
  //You may need to adjust DIST_MULTIPLIER especialy with non zero Julia seed
//  return (DIST_MULTIPLIER*de_PZshape(p-PK_Offset)/abs(DEfactor)-PK_DEoffset);
  return (DIST_MULTIPLIER*de_KIFSMenger(p-PK_Offset)/abs(DEfactor)-PK_DEoffset);
}
DECLARE_DE(de_PKlein)

// Compute the color.
vec3 c_PKlein(vec3 p) {
   //Just scale=1 Julia box
  float r2=dot(p,p);
  //float DEfactor=1.;
  vec3  col=vec3(0.0);
  float rmin=10000.0;

  for(int i=0;i<color_iters;i++){
    vec3 p1=clamp(p, -PK_BSize, PK_BSize)*2. - p;  // Box folding
    col+=abs(p-p1);
    p=p1;
    r2=dot(p,p);  // Inversion
    float k=max(PK_CSize/r2,1.);
    p*=k; //DEfactor*=k;
    p+=PK_CVector;  // julia seed
    r2=dot(p,p);
    rmin=min(rmin,r2);
  }
  return mix(vec3(sqrt(rmin)),(vec3(0.5)+sin(PK_Color*col.z)*0.5),PK_CMix);
}
DECLARE_COLORING(c_PKlein)

#define rotationVector par[3]
#define rotationAngle par[4].x  // { min=-5 max=5 step=.01}

float de_mandelbox(vec3 pos) {
  float minRad2 = clamp(float(MB_MINRAD2), float(1.0e-9), float(1.0));
  vec4 scale = vec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  
  float s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
  for (int i=0; i<iters; i++) s*= ds;
  float absScalePowIters = s;
  
  float csat = cos(rotationAngle);
  float ssat = sin(rotationAngle);
  float usat = 1.0 - csat;
  vec3 u = normalize(rotationVector);
  mat3 rotationMatrix = mat3(
    u.x*u.x*usat + csat,     u.x*u.y*usat - u.z*ssat, u.x*u.z*usat + u.y*ssat,
    u.y*u.x*usat + u.z*ssat, u.y*u.y*usat + csat,     u.y*u.z*usat - u.x*ssat,
    u.z*u.x*usat - u.y*ssat, u.z*u.y*usat + u.x*ssat, u.z*u.z*usat + csat
    );

  vec4 p = vec4(pos, 1.0);  // p.w is the distance estimate
  vec4 p0;
  if (julia) p0 = vec4(JuliaVector, 1.0); else p0 = p;
  for (int i=0; i<iters; i++) {
    p = vec4(rotationMatrix * p.xyz, p.w);
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
    if (r2 > 100.0) break;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
     - absScalePowIters) * DIST_MULTIPLIER;
//    return (de_KIFSMenger(p.xyz)/abs(p.w))*DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

#ifdef _FAKE_GLSL_

double de_mandelbox_64(dvec3 pos) {
  double minRad2 = clamp(double(MB_MINRAD2), 1.0e-9, 1.0);
  dvec4 scale = dvec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  
  double s = abs(MB_SCALE), ds = 1.0 / abs(MB_SCALE);
  for (int i=0; i<iters; i++) s*= ds;
  double absScalePowIters = s;
  
  float csat = cos(-rotationAngle);  // Note negate: messed up host vs. shader orientation?!
  float ssat = sin(-rotationAngle);
  float usat = 1.0f - csat;
  vec3 u = normalize(rotationVector);
  mat3 rotationMatrix = mat3(
    u.x*u.x*usat + csat,     u.x*u.y*usat - u.z*ssat, u.x*u.z*usat + u.y*ssat,
    u.y*u.x*usat + u.z*ssat, u.y*u.y*usat + csat,     u.y*u.z*usat - u.x*ssat,
    u.z*u.x*usat - u.y*ssat, u.z*u.y*usat + u.x*ssat, u.z*u.z*usat + csat
    );

  dvec4 p = dvec4(pos, 1.0);  // p.w is the distance estimate
  dvec4 p0;
  if (julia)
    p0 = dvec4(JuliaVector.x, JuliaVector.y, JuliaVector.z, 1.0);
  else
    p0 = p;
  for (int i=0; i<iters; i++) {
    p = dvec4(rotationMatrix * dvec3(p.xyz), p.w);
    p = dvec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    double r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;

	//if (r2 > 100.0) break;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
			- absScalePowIters) * 0.95 * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox_64)

double de_Z0_64(dvec3 pos) {
  return abs(pos.z * .3);
}
DECLARE_DE(de_Z0_64)

#endif  // _FAKE_GLSL_

// Infinite construction menger sphere sponge.
float de_ssponge(vec3 pos) {
  float k = SS_INITIAL_K;
  float d = -100.0;
  for(int i=0; i<iters; i++) {
    vec3 x = mod(pos * k, SS_MOD) - .5 * SS_MOD;
    float r = length(x);
    float d1 = (SS_INTRA_D - r) / k;
    d = max(d, d1);
    k *= SS_SCALE;
  }
  return d;
}
DECLARE_DE(de_ssponge)

// DE for cylinder
float de_cylinder(vec3 pos) {
  return length(pos.xz - CYL_Vector.xy) - CYL_Vector.z;
}
DECLARE_DE(de_cylinder)

// DE for sphere
float de_sphere(vec3 pos) {
  return length(pos) - SPHERE_SZ;
}
DECLARE_DE(de_sphere)

float de_combi(vec3 z0) {
  // Use min() for union, max() for intersection of shapes.
  // max(-a, b) for subtraction.
  //return min(max(de_ssponge(z0),de_mandelbox(z0)), de_menger(z0));
  //return max(-de_cylinder(z0), de_mandelbox(z0));
  return min(de_sphere(z0), de_mandelbox(z0));
}
DECLARE_DE(de_combi)

// Compute the color at `pos`.
vec3 c_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, float(1.0e-9), float(1.0));
  vec3 scale = vec3(MB_SCALE, MB_SCALE, MB_SCALE) / minRad2;
  vec3 p = pos;
  vec3 p0;
  float trap = 1.0;

  if (julia) p0 = JuliaVector; else p0 = p;
  for (int i=0; i<color_iters; i++) {
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    float r2 = dot(p, p);
    p *= clamp(float(max(minRad2/r2, minRad2)), float(0.0), float(1.0));
    p = p*scale + p0;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);
  // return texture2DLod(bg_texture, c, 0.0).xyz;
  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}
DECLARE_COLORING(c_mandelbox)

// --- Refectoids DE, imported here so boxplorer.cc can haz DE; more fun to fly around with.
// --- note in refectoids.cfg, add d d_PKlein to let .exe know which DE to use.
// Compute the distance from `pos` to the PKlein basic shape.
float d_PZshape(vec3 p) {
  float rxy = sign(par[9].y)*(length(p.xy)-abs(par[9].y));

  float TThickness = par[5].x;
  float Zmult = par[8].y;
  float Ziter = par[8].x;
 
  for(int i=0; i<int(Ziter); i++) p.z=2.*clamp(p.z, -Zmult, Zmult)-p.z;

// This abs() creates the nice holes but also causes banding.
// For movement-only DE, holes are just what we want ;-)
  return max(rxy,abs(length(p.xy)*p.z-TThickness) / sqrt(dot(p,p)+abs(TThickness)));
 //  return max(rxy,   -(length(p.xy)*p.z-TThickness) / sqrt(dot(p,p)+abs(TThickness)));
}

// Compute the distance from `pos` to the PKlein.
float d_PKlein(vec3 p) {
   //Just scale=1 Julia box
  float r2=dot(p,p);
  float DEfactor=1.;

  float rDIST_MULTIPLIER = par[9].x;
  vec3 CSize = vec3(par[1].y,par[1].x,par[2].y);
  float Size = par[0].y;
  vec3 C = vec3(par[2].x,par[3].y,par[3].x);
  vec3 Offset = vec3(par[4].y,par[4].x,par[5].y);
  float DEoffset = par[0].x;

	for(int i=0;i<iters;i++){
		//Box folding
		p=clamp(p, -CSize, CSize) * 2.0 -p;
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
	return (rDIST_MULTIPLIER*d_PZshape(p-Offset)/abs(DEfactor)-DEoffset);
}
DECLARE_DE(d_PKlein)
// -- end of Reflectoids DE

const float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
//  vec2 Eps = vec2(0, max(d_pos, normal_eps));
//  vec2 Eps = vec2(0, d_pos);
  vec2 Eps = vec2(0, max(normal_eps, d_pos));
  return normalize(vec3(
  // 3-tap central differences, error = O(eps^2)
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)
  ));
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor, vec3 specular) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(float(max( float(dot(normal, halfLV)), 0.0 )), float(32.0));
  float dif = dot(normal, light) * 0.5 + 0.75;
  return diffuseColor*dif + specular*spe;
}

// Mix reflected ray contribution.
vec3 mix_reflection(vec3 normal, vec3 view, vec3 baseColor, vec3 reflectionColor, float factor) {
  // An attempt at Oren-Nayar reflectance
  float alpha = acos(dot(normal, -view));
  float s2 = SIGMA * SIGMA;
  float A = 1. - .5 * (s2 / (s2 + .33));
  float B = .45 * s2 / (s2 + .09);
  float rho = ALBEDO;
  float li = 1.0;  // incident intensity, pegged at 1 for now.
  float a = max(alpha, alpha);  // incident and reflection angles are same..
  float b = min(alpha, alpha);
  float ri = rho * cos(alpha) * ( A + (B * max(0., cos(alpha - alpha)) * sin(a) * tan(b))) * li;
  return mix(baseColor, reflectionColor, abs(factor * ri));
}

// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, float totalD, float m_dist, float side) {
  float ao_ed = ao_eps; //totalD*ao_eps/m_dist;
  float ao = 1.0, w = ao_strength/ao_ed;
  float dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-abs(D)) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(ao, float(0.0), float(1.0));
}

#ifndef _FAKE_GLSL_
// Intersect direction ray w/ large encapsulating sphere.
// Sphere map texture onto it.
vec3 background_color(in vec3 vp) {
  const vec3 vn = vec3(0.0, 1.0, 0.0);
  const vec3 ve = vec3(1.0, 0.0, 0.0);
  float phi = acos(-dot(vn, vp));
  float v = phi / PI;
  float theta = (acos(dot(vp, ve) / sin(phi))) / (2. * PI);
  float u;
  if (dot(cross(vn, ve), vp) > 0.) {
    u = theta;
  } else {
    u = 1. - theta;
  }
  //return texture2DLod(bg_texture, vec2(u+time/10.0,v+time/10.0), BG_BLUR).xyz;
  return texture2DLod(iChannel0, vec2(u+time/40.0,v), BG_BLUR).xyz * 2.0;
  //return texture2D(bg_texture, vec2(u+time/10.0,v)).xyz;
}
#else
// TODO: add texture2DLod to fake glsl
vec3 background_color(vec3 vp) { return backgroundColor; }
#endif

// Intersect the view ray with the fractal using raymarching.
// returns # steps
int rayMarch(vec3 p, vec3 dp, INOUT(float,totalD), float side, INOUT(float,m_dist), float m_zoom) {
  int steps;
  float D = 0.0;
  for (steps = 0; steps < max_steps; ++steps) {
   // float D = (side * d(p + dp * totalD) - totalD * m_dist) / (1.0 + m_dist);
   D = side * d(p + dp * totalD);
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist = m_zoom * totalD;
  }
#if 1
 if (D < m_dist)
for (int i = 0; i<5; ++i) {
  totalD += D - m_dist;
  m_dist = m_zoom * totalD;
  D = abs(side * d(p + dp * totalD));
}
#endif
  return steps;
}

// Trace from point p back to light(s). Return received direct light.
vec3 shade(vec3 p, vec3 l1, float side, float m_zoom) {
  vec3 dir = l1 - p;
  float l1_dist = length(dir);
  dir = normalize(dir);
  float fudge = 2.0 * min_dist;
  p += dir * fudge;
  float totalD = 0.0;
  float m_dist = max(min_dist, m_zoom * totalD);  // reset m_dist, fresh ray
  rayMarch(p, dir, totalD, side, m_dist, m_zoom);  // trace towards light
  if (totalD >= l1_dist - fudge) {
    float falloff = pow(1.0 + l1_dist, 2.0);
	return vec3(0.,0.,1.0) * clamp(2.0 / falloff, 0.0, 1.0);
  }
  return vec3(0.0, 0.0, 0.0);
}

// return color and ratio to mix in pure light from L1_Vector (direct visibility).
vec4 lightBulb(vec3 x2, vec3 dp, float totalD) {
  vec3 x1 = x2 - dp * totalD;
  vec3 x0 = L1_Vector;
  float t = -dot(x1 - x0, x2 - x1) / dot(x2 - x1, x2 - x1);
  if (t < 0.0 || t > 1.0) return vec4(0.0);  // not near this segment.
  float d = length(cross(x0 - x1, x0 - x2)) / length(x2 - x1);
  float delta = L1_Size - d;
  if (delta < 0.0) return vec4(0.0);  // larger than light radius
  delta /= L1_Size;
  return vec4(
    clamp(1.3*delta, 0.0, 1.0),
    clamp(1.3*delta, 0.0, 1.0),
    1.0,
    clamp(3.0*delta, 0.0, 0.99 /*1.0 causes artifacts.. dunno why*/));
}

// Get base color at p, plus Blinn_Phong and ambient occulusion.
vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist, float side, float m_zoom) {
  vec3 col = c(p);
  col = blinn_phong(n, -dp, normalize(vec3(1.0,.6,0.7)+dp), col,
  //mix(specularColor, col, pow(abs((m_dist-min_dist))/(m_dist+min_dist), float(2.5))));
     specularColor);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist, side));
  col += shade(p, L1_Vector, side, m_zoom);
  return col;
}


#if 0
// 3d noise
vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}
vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}
vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}
vec4 taylorInvSqrt(vec4 r) {
  return 1.79284291400159 - 0.85373472095314 * r;
}
float snoise(vec3 v) {
  const vec2 C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);
  vec3 i = floor(v + dot(v, C.yyy) );
  vec3 x0 = v - i + dot(i, C.xxx) ;
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy;
  vec3 x3 = x0 - D.yyy;
  i = mod289(i);
  vec4 p = permute( permute( permute(
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 ))
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));
  float n_ = 0.142857142857;
  vec3 ns = n_ * D.wyz - D.xzx;
  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);
  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );
  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);
  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));
  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;
  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1),
                                dot(p2,x2), dot(p3,x3) ) );
}
#endif

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){ return fract(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547); }

void main() {
  vec3 eye_out, dp; 

  if (!setup_ray(eye, dir, eye_out, dp)) return;

  float m_zoom = zoom * .5 / xres;  // screen error at dist 1.
  float noise = pnoise(gl_FragCoord.xy);

  vec3 p = eye_out;
  float totalD = d(p);
  float side = sign(totalD);
  totalD *= side * noise;  // Randomize first step.
  float m_dist = m_zoom * totalD;

  int steps = 0;  // number of marching steps we've taken for this ray.
  float colFactor = SHINE;

  // March first ray.
  steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);
  vec3 rayCol, n;
  vec4 light = lightBulb(p + dp * totalD, dp, totalD);
  if (totalD < MAX_DIST) {
    p += dp * totalD;
    n = normal(p, m_dist * .5);
    rayCol = rayColor(p, dp, n, totalD, m_dist, side, m_zoom);
    rayCol = mix(rayCol, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);
  } else {
    rayCol = background_color(dp);
    //totalD = 0.;
  }

  vec3 finalCol = rayCol;
  float firstD = totalD;

if (firstD != 0.0)
// uncomment to have only sphere be reflective.
//  if (de_sphere(p) < de_mandelbox(p))
{
//     n = normalize(p);

  // March reflected ray a couple of times.
  for (int ray = 1; ray < nrays &&
                    totalD < MAX_DIST &&
                    colFactor > 0.0; ++ray) {
    dp = reflect(dp, n);  // reflect view direction
    p += dp * (-totalD + (3.+noise) * m_dist);  // reproject eye

    float oldTotalD = totalD;  // only trace lightbulbs on actual reflected path..
    steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);

    vec4 rlight = lightBulb(p + dp * totalD, dp, totalD - oldTotalD);
    if (totalD < MAX_DIST) {
      p += dp * totalD;
      n = normal(p, m_dist * .5);
      rayCol = rayColor(p, dp, n, totalD, m_dist, side, m_zoom);
      rayCol = mix(rayCol, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);
    } else {
      rayCol = background_color(dp);
    }
    finalCol = mix_reflection(n, -dp, finalCol, rayCol, colFactor);
    finalCol = mix(finalCol, rlight.xyz, rlight.w * colFactor);
    colFactor *= (SHINE * SHINE);  // reflection drop-off.
  }
}

  // draw lights, if any on primary ray.
  finalCol = mix(finalCol, light.xyz, light.w);

  // gamma
  //finalCol = pow(clamp(finalCol, 0.0, 1.0), vec3(0.45));

  write_pixel(dir, firstD, finalCol);
}
