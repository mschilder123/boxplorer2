// http://www.fractalforums.com/general-discussion-b77/solids-many-many-solids/msg43794/#msg43794
// from .frag by Knighty

#include "setup.inc"
#line 7
 
#define INOUT(a,b) inout a b
#define MAX_DIST 10.0
#define ULP 0.000000059604644775390625
#ifndef PI
#define PI 3.14159265
#endif

uniform vec3 par[20];
uniform int nrays;        // {min=1 max=10} # of ray bounces.

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters,    // Number of fractal iterations.
  color_iters,        // Number of fractal iterations for coloring.
  max_steps;          // Maximum raymarching steps.

// Reflection params
#define SHINE par[9].x  // {min=0 max=1 step=.01} reflection contribution drop-off
#define SIGMA par[9].y  //{min=0 max=1 step=.01} surface roughness
#define ALBEDO par[9].z  //{min=0 max=1 step=.01} surface reflection

// Background blur
#define BG_BLUR par[1].z  // {min=0 max=8 step=.1}

//#group Hyperbolic-tesselation

// UVWT 'barycentric' coordinate for the 'principal' node
#define U par[1].x  // {min=0 max=1 step=.01}
#define V par[1].y  // {min=0 max=1 step=.01}
#define W par[2].x  // {min=0 max=1 step=.01}
#define T par[2].y  // {min=0 max=1 step=.01}

//vertex radius
#define VRadius par[3].x  // {min=0 max=1 step=.01}

//segments radius
#define SRadius par[3].y  // {min=0 max=1 step=.01}

//cutting sphere radius
#define CSphRad par[4].y  // {min=0 max=1 step=.01}

#define RotAngle par[4].x
#define RotVector par[14]

//#group HTess-Color
#define segAColor par[19]
#define segBColor par[18]
#define segCColor par[17]
#define segDColor par[16]
#define verticesColor par[15]

#define backgroundColor par[10]
#define specularColor par[11]
#define glowColor par[12]
#define aoColor par[13]

#define L1_Vector par[8]
#define L1_Size par[1].z  // {min=0 max=1 step=.01}

vec4 nc,nd,p;
float cVR,sVR,cSR,sSR,cRA,sRA;
float hdot(vec4 a, vec4 b){//dot product for Minkowski space.
  return dot(a.xyz,b.xyz)-a.w*b.w;
}

vec4 hnormalizew(vec4 v){//normalization of (timelike) vectors in Minkowski space.
  float l=1./sqrt(v.w*v.w-dot(v.xyz,v.xyz));
  return v*l;
}

float cosh(float val) {
  float tmp = exp(val);
  float cosH = (tmp + 1.0 / tmp) / 2.0;
  return cosH;
}

float sinh(float val) {
  float tmp = exp(val);
  float sinH = (tmp - 1.0 / tmp) / 2.0;
  return sinH;
}

void init() {
  float cospin=cos(PI/5.);
  float scospin=sqrt(4./3.*cospin*cospin-3./4.);

  //na and nb are simply vec4(1.,0.,0.,0.) and vec4(0.,1.,0.,0.) respectively
  nc=0.5*vec4(0,-1,sqrt(3.),0.);
  nd=vec4(-0.5,-cospin,-cospin/sqrt(3.),-scospin);

  vec4 pabc,pbdc,pcda,pdba;
  pabc=vec4(0.,0.,0.,0.5*sqrt(3.));
  pbdc=0.5*sqrt(3.)*vec4(scospin,0.,0.,0.5);
  pcda=vec4(0.,0.5*sqrt(3.)*scospin,0.5*scospin,cospin*2./sqrt(3.));
  pdba=vec4(0.,0.,scospin,cospin/sqrt(3.));

  p=hnormalizew(U*pabc+V*pbdc+W*pcda+T*pdba);

  cVR=cosh(VRadius);sVR=sinh(VRadius);
  cSR=cosh(SRadius);sSR=sinh(SRadius);
  cRA=cosh(RotAngle);sRA=-sinh(RotAngle);
}

vec4 Rotate(vec4 p){
  //this is a (hyperbolic) rotation (that is, a boost) on the plane defined by RotVector and w axis
  //We do not need more because the remaining 3 rotation are in our 3D space
  //That would be redundant.
  //This rotation is equivalent to a translation inside the hyperbolic space when the camera is at 0,0,0
  vec4 p1=p;
  vec3 rv;
  rv=normalize(RotVector);
  float vp=dot(rv,p.xyz);
  p1.xyz+=rv*(vp*(cRA-1.)+p.w*sRA);
  p1.w+=vp*sRA+p.w*(cRA-1.);
  return p1;
}

vec4 fold(vec4 pos) {//beside using minkowski dot product, its exactly the same as for euclidean space
  for(int i=0;i<iters;i++){
    pos.xy=abs(pos.xy);
    float t=-2.*min(0.,hdot(pos,nc)); pos+=t*nc;
    t=-2.*min(0.,hdot(pos,nd)); pos+=t*nd;
  }
  return pos;
}

float DD(float ca, float sa, float r){//converts hyperbolic distance to distance in projection flat space. ca and sa are the hyperbolic cosine and sine of the hyperbolic distance which is an "angle".
  return (2.*r*ca+(1.+r*r)*sa)/((1.+r*r)*ca+2.*r*sa+1.-r*r)-r;
}

float dist2Vertex(vec4 z, float r){
  float ca=-hdot(z,p), sa=0.5*sqrt(-hdot(p-z,p-z)*hdot(p+z,p+z));

  return DD(ca*cVR-sa*sVR,sa*cVR-ca*sVR,r);
}

float dist2Segment(vec4 z, vec4 n, float r){
  //pmin is the orthogonal projection of z onto the plane defined by p and n
  //then pmin is projected onto the unit sphere
  float zn=hdot(z,n),zp=hdot(z,p),np=hdot(n,p);
  float det=-1./(1.+np*np);
  float alpha=det*(zp-zn*np), beta=det*(-np*zp-zn);
  vec4 pmin=hnormalizew(alpha*p+min(0.,beta)*n);
  //ca and sa are the hyperbolic cosine and sine of the "angle" between z and pmin. This is the distance in hyperbolic space.
  float ca=-hdot(z,pmin), sa=0.5*sqrt(-hdot(pmin-z,pmin-z)*hdot(pmin+z,pmin+z));
  return DD(ca*cSR-sa*sSR,sa*cSR-ca*sSR,r);//we subtract the width of the sgment before conversion
}
//it is possible to compute the distance to a face just as for segments: pmin will be the orthogonal projection
// of z onto the 3-plane defined by p and two n's (na and nb, na and nc, na and and, nb and nd... and so on).
//that involves solving a system of 3 linear equations.
//it's not implemented here because it is better with transparency

float dist2Segments(vec4 z, float r){
  float da=dist2Segment(z, vec4(1.,0.,0.,0.), r);
  float db=dist2Segment(z, vec4(0.,1.,0.,0.), r);
  float dc=dist2Segment(z, nc, r);
  float dd=dist2Segment(z, nd, r);

  return min(min(da,db),min(dc,dd));
}

float d_verts(vec3 pos) {
  float r=length(pos);
  vec4 z4=vec4(2.*pos,1.+r*r)*1./(1.-r*r);//Inverse stereographic projection of pos: z4 lies onto the unit 3-parabolid of revolution around w axis centered at 0.
  z4=Rotate(z4);
  z4=fold(z4);
  return min(dist2Vertex(z4,r),dist2Segments(z4,r)) *.95;
}

float d_sphere(vec3 pos) {
  float r=length(pos);
  return r-CSphRad;
}

float d(vec3 pos) {
  float ds = d_sphere(pos);
  if (ds > 0.0) return ds;
  return min(-ds, d_verts(pos));
}

vec3 c(vec3 pos){
  float ds = d_sphere(pos);
  if (ds > 0.0 || -ds < d_verts(pos)) return vec3(.7,.7,.7);
  float r=length(pos);
  vec4 z4=vec4(2.*pos,1.+r*r)*1./(1.-r*r);
  z4=Rotate(z4);
  z4=fold(z4);
  float da=dist2Segment(z4, vec4(1.,0.,0.,0.), r);
  float db=dist2Segment(z4, vec4(0.,1.,0.,0.), r);
  float dc=dist2Segment(z4, nc, r);
  float dd=dist2Segment(z4, nd, r);
  float dv=dist2Vertex(z4,r);
  float d=min(min(min(da,db),min(dc,dd)),dv);
  vec3 color=segAColor;
  if(d==db) color=segBColor;
  if(d==dc) color=segCColor;
  if(d==dd) color=segDColor;
  if(d==dv) color=verticesColor;
  return color;
}

const float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
//  vec2 Eps = vec2(0, max(d_pos, normal_eps));
  vec2 Eps = vec2(0, max(d_pos,normal_eps));
//  vec2 Eps = vec2(0, normal_eps);
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
  float ao_ed = m_dist; //ao_eps; //totalD*ao_eps/m_dist;
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
uniform sampler2D bg_texture;
uniform int use_bg_texture;
vec3 background_color(in vec3 vp) {
  if (use_bg_texture == 0) return backgroundColor;
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
  //return texture2DLod(bg_texture, vec2(u+time/10.0,v), BG_BLUR).xyz;
  return texture2D(bg_texture, vec2(u+time/10.0,v)).xyz;
}
#else
// TODO: add texture2DLod to fake glsl
vec3 background_color(vec3 vp) { return backgroundColor; }
#endif

// Intersect the view ray with the fractal using raymarching.
// returns # steps
int rayMarch(vec3 p, vec3 dp, INOUT(float,totalD), float side, INOUT(float,m_dist), float m_zoom) {
  int steps;
  for (steps = 0; steps < max_steps; ++steps) {
    // float D = (side * d(p + dp * totalD) - totalD * m_dist) / (1.0 + m_dist);
    float D = side * d(p + dp * totalD);
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist = m_zoom * totalD;
  }
  return steps;
}

// Trace from point p back to light(s). Return received direct light.
vec3 shade(vec3 p, vec3 l1, float side, float m_zoom, vec3 n, vec3 col, vec3 dp) {
  if (L1_Size == 0.0) return col;
  vec3 dir = l1 - p;
  float l1_dist = length(dir);
  dir = normalize(dir);
  float fudge = 2.0 * min_dist;
  p += dir * fudge;
  float totalD = 0.0;
  float m_dist = max(min_dist, m_zoom * totalD);  // reset m_dist, fresh ray
  rayMarch(p, dir, totalD, side, m_dist, m_zoom);  // trace towards light

  float falloff = pow(1.0 + l1_dist, 3.0);
  vec3 c = vec3(0.3,0.3,0.4) * clamp(2.0 / falloff, 0.0, 1.0);
  c = blinn_phong(n, -dp, dir, col, c);
  if (totalD < l1_dist + fudge) c = c * .5;
  return c;
}

vec3 getLight1() {
  vec3 l = L1_Vector;
  return l * vec3(0,cos(time),sin(time));  // circle around
}

// return color and ratio to mix in pure light from L1_Vector (direct visibility).
vec4 lightBulb(vec3 x2, vec3 dp, float totalD) {
  vec3 x1 = x2 - dp * totalD;
  vec3 x0 = getLight1();
  float t = -dot(x1 - x0, x2 - x1) / dot(x2 - x1, x2 - x1);
  if (t <= 0.0 || t >= 1.0) return vec4(0.0);  // not near this segment.
  float d = length(cross(x0 - x1, x0 - x2)) / length(x2 - x1);
  if (d > L1_Size) return vec4(0.0);  // larger than light radius
  return vec4(
    clamp(1.3*(L1_Size-d)/L1_Size, 0.0, 1.0),
    clamp(1.3*(L1_Size-d)/L1_Size, 0.0, 1.0),
    clamp(1.3*(L1_Size-d)/L1_Size, 0.0, 1.0),
    clamp(pow(1.5*(L1_Size-d)/L1_Size, 3.0), 0.0, 0.99));
}

// Get base color at p, plus Blinn_phing and ambient occulusion.
vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist, float side, float m_zoom) {
  vec3 col = c(p);
  col = blinn_phong(n, -dp, normalize(vec3(1.0,.6,0.7)+dp), col, specularColor);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist, side));
  col = shade(p, getLight1(), side, m_zoom, n, col, dp);
  return col;
}

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){ return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

void main() {
  vec3 eye_in, dp;
  if (!setup_ray(eye, dir ,eye_in, dp)) return;

  init();

  float m_zoom = zoom * .5 / xres;  // screen error at dist 1.
  float noise = pnoise(gl_FragCoord.xy);

  vec3 p = eye_in;
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
  bool sphereHit = false;
  if (totalD < MAX_DIST) {
    p += dp * totalD;
    float ds = d_sphere(p);
    if (ds > 0.0 || -ds < d_verts(p)) {
      sphereHit = true;
      n = normalize(p);  // easy accurate normal for sphere hit
    } else {
      n = normal(p, m_dist * .5);
    }
    rayCol = rayColor(p, dp, n, totalD, m_dist, side, m_zoom);
    rayCol = mix(rayCol, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);
  } else {
    rayCol = background_color(dp);
  }

  float firstD = totalD;
  vec3 finalCol = rayCol;

  // March reflected ray a couple of times.
  for (int ray = 1; ray < nrays &&
          sphereHit &&
          totalD < MAX_DIST &&
          colFactor > 0.0; ++ray) {
    sphereHit = false;
    dp = reflect(dp, n);  // reflect view direction
    p += dp * (-totalD + (1.5+noise) * m_dist);  // reproject eye

    float oldTotalD = totalD;  // only trace lightbulbs on actual reflected path..
    steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);

    vec4 rlight = lightBulb(p + dp * totalD, dp, totalD - oldTotalD);
    if (totalD < MAX_DIST) {
      p += dp * totalD;
      float ds = d_sphere(p);
      if (ds > 0.0 || -ds < d_verts(p)) {
        sphereHit = true;
        n = normalize(p);  // easy accurate normal for sphere hit
      } else {
        n = normal(p, m_dist * .5);
      }
      rayCol = rayColor(p, dp, n, totalD, m_dist, side, m_zoom);
      rayCol = mix(rayCol, glowColor, (float(steps)+noise)/float(max_steps) * glow_strength);
    } else {
      rayCol = background_color(dp);
    }
    finalCol = mix_reflection(n, -dp, finalCol, rayCol, colFactor);
    finalCol = mix(finalCol, rlight.xyz, rlight.w * colFactor);
    colFactor *= SHINE * SHINE;  // reflection drop-off.
  }

  // draw lights, if any on primary ray.
  finalCol = mix(finalCol, light.xyz, light.w);

  write_pixel(dir, firstD, finalCol);
}
