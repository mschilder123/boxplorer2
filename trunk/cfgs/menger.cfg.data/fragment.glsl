// menger and then some shader.
// Original shader by rrrola for mandelbox
// bermarte, marius: formulae from Knighty
// marius: refactored w/ reflections, background dome, ssponge, combi etc.
// marius: massaged so can be compiled as plain C++.

#ifndef _FAKE_GLSL_

#define DECLARE_DE(x)
#define DECLARE_COLORING(x)
#define INOUT(a,b) inout a b

// distance estimator func
#ifndef d
#define d de_menger // PKlein,combi,menger,mandelbox,ssponge
#endif

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

// Camera position and direction.
varying vec3 eye, dir;
varying float zoom;
uniform float xres;
uniform float time;

// Interactive parameters.
uniform vec3 par[20];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform float speed;  // eye separation really.

uniform int iters;        // {min=1 max=1000} Number of fractal iterations.
uniform int color_iters;  // {min=0 max=1000} Coloration iterations.
uniform int max_steps;    // {min=1 max=1000} Maximum raymarching steps.
uniform int nrays;        // {min=1 max=10} # of ray bounces.

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
#define ME_SIZE par[1].y // {min=0 max=5 step=.01} size of menger box

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
#define PK_CMix par[7].x
#define PK_DEoffset par[0].x

#define CYL_Vector par[17]

#define L1_Vector par[19]
#define L1_Size par[18].x  // {min=0.001 max=1.0 step=.001}

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
  vec3 ap=p+vec3(1.);

  for(int i=0;i<iters;i++) {
    ap=p;
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
  return (DIST_MULTIPLIER*de_PZshape(p-PK_Offset)/abs(DEfactor)-PK_DEoffset);
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
DECLARE_DE(de_menger)

vec3 c_menger(vec3 p) {
  return surfaceColor1;  // Boring but reflections make it interesting.
}
DECLARE_COLORING(c_menger)

float de_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  vec4 scale = vec4(MB_SCALE, MB_SCALE, MB_SCALE, abs(MB_SCALE)) / minRad2;
  vec4 p = vec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  for (int i=0; i<iters; i++) {
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale + p0;
  }
  return ((length(p.xyz) - abs(MB_SCALE - 1.0)) / p.w
           - pow(abs(MB_SCALE), float(1-iters))) * DIST_MULTIPLIER;
}
DECLARE_DE(de_mandelbox)

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

float de_combi(vec3 z0) {
  // Use min() for union, max() for intersection of shapes.
  // max(-a, b) for subtraction.
  //return min(max(de_ssponge(z0),de_mandelbox(z0)), de_menger(z0));
  return max(-de_cylinder(z0), de_mandelbox(z0));
}
DECLARE_DE(de_combi)

const float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec2 Eps = vec2(0, max(d_pos, normal_eps));
//  vec2 Eps = vec2(0, d_pos);
  return normalize(vec3(
  // 3-tap central differences, error = O(eps^2)
    -d(pos-Eps.yxx)+d(pos+Eps.yxx),
    -d(pos-Eps.xyx)+d(pos+Eps.xyx),
    -d(pos-Eps.xxy)+d(pos+Eps.xxy)
  ));
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), float(32.0));
  float dif = dot(normal, light) * 0.5 + 0.75;
  return diffuseColor*dif + specularColor*spe;
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
  float ao_ed = totalD*ao_eps/m_dist;
  float ao = 1.0, w = ao_strength/ao_ed;
  float dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
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
  return texture2DLod(bg_texture, vec2(u+time/10.0,v), BG_BLUR).xyz;
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
    float D = (side * d(p + dp * totalD) - totalD * m_dist) / (1.0 + m_dist);
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist = max(min_dist, m_zoom * totalD);
  }
  return steps;
}

// Compute the color at `pos`.
vec3 c_mandelbox(vec3 pos) {
  float minRad2 = clamp(MB_MINRAD2, 1.0e-9, 1.0);
  vec3 scale = vec3(MB_SCALE, MB_SCALE, MB_SCALE) / minRad2;
  vec3 p = pos, p0 = p;
  float trap = 1.0;

  for (int i=0; i<color_iters; i++) {
    p = clamp(p, -1.0, 1.0) * 2.0 - p;
    float r2 = dot(p, p);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
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
  if (d > L1_Size) return vec4(0.0);  // larger than light radius
  return vec4(
    clamp(1.3*(L1_Size-d)/L1_Size, 0.0, 1.0),
    clamp(1.3*(L1_Size-d)/L1_Size, 0.0, 1.0),
    1.0,
	clamp(3.0*(L1_Size-d)/L1_Size, 0.0, 1.0));
}

// Get base color at p, plus Blinn_phing and ambient occulusion.
vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist, float side, float m_zoom) {
  vec3 col = c(p);
  col = blinn_phong(n, -dp, normalize(vec3(1.0,.6,0.7)+dp), col);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist, side));
  col += shade(p, L1_Vector, side, m_zoom);
  return col;
}

void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;

#ifndef _FAKE_GLSL_
  // TODO: add gl_FragCoord.xy and gl_ModelViewMatrix to fake glsl
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);
#endif

  vec3 p = eye_in, dp = normalize(dir);
  float m_zoom = /*length(dir) * */ zoom * .25 / xres;  // screen error at dist 1.

  float totalD = d(p);
  float side = sign(totalD);
  totalD *= side * .5;  // start with conservative step.
  float m_dist = max(min_dist, m_zoom * totalD);

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
    rayCol = mix(rayCol, glowColor, float(steps)/float(max_steps) * glow_strength);
  } else {
    rayCol = background_color(dp);
  }

  float firstD = totalD;
  vec3 finalCol = rayCol;

  // March reflected ray a couple of times.
  for (int ray = 1; ray < nrays &&
                    totalD < MAX_DIST &&
                    colFactor > 0.0; ++ray) {
    dp = reflect(dp, n);  // reflect view direction
    p += dp * (-totalD + 2.0 * m_dist);  // reproject eye

	float oldTotalD = totalD;  // only trace lightbulbs on actual reflected path..
    steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);

	vec4 rlight = lightBulb(p + dp * totalD, dp, totalD - oldTotalD);
    if (totalD < MAX_DIST) {
      p += dp * totalD;
      n = normal(p, m_dist * .5);
      rayCol = rayColor(p, dp, n, totalD, m_dist, side, m_zoom);
      rayCol = mix(rayCol, glowColor, float(steps)/float(max_steps) * glow_strength);
    } else {
      rayCol = background_color(dp);
    }
    finalCol = mix_reflection(n, -dp, finalCol, rayCol, colFactor);
	finalCol = mix(finalCol, rlight.xyz, rlight.w * colFactor);
    colFactor *= SHINE * SHINE;  // reflection drop-off.
  }

  // draw lights, if any on primary ray.
  finalCol = mix(finalCol, light.xyz, light.w);

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(firstD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(finalCol, depth);
}

