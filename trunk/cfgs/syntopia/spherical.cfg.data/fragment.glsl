// http://blog.hvidtfeldts.net/index.php/2012/03/spherical-worlds/
// Adapted for boxplorer2 by marius.

#define INOUT(a,b) inout a b
#define MAX_DIST 20.0
#define ULP 0.000000059604644775390625
#ifndef PI
#define PI 3.14159265
#endif

uniform vec3 par[20];
uniform float speed;
varying float zoom;
uniform float xres;
uniform float time;
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

#define backgroundColor par[10]
#define specularColor par[11]
#define glowColor par[12]
#define aoColor par[13]

#define L1_Vector par[8]
#define L1_Size par[1].z  // {min=0 max=1 step=.01}

varying vec3 eye, dir;

#if __VERSION__ < 130
float cosh(float val)
{
  float tmp = exp(val);
  float cosH = (tmp + 1.0 / tmp) / 2.0;
  return cosH;
}

float tanh(float val)
{
  float tmp = exp(val);
  float tanH = (tmp - 1.0 / tmp) / (tmp + 1.0 / tmp);
  return tanH;
}

float sinh(float val)
{
  float tmp = exp(val);
  float sinH = (tmp - 1.0 / tmp) / 2.0;
  return sinH;
}
#endif

vec2 cMul(vec2 a, vec2 b) {
	return vec2( a.x*b.x -  a.y*b.y,a.x*b.y + a.y * b.x);
}

vec2 cPower(vec2 z, float n) {
	float r2 = dot(z,z);
	return pow(r2,n/2.0)*vec2(cos(n*atan(z.y/z.x)),sin(n*atan(z.y/z.x)));
}

vec2 cInverse(vec2 a) {
  return vec2(a.x,-a.y)/dot(a,a);
}

vec2 cDiv(vec2 a, vec2 b) {
	return cMul(a,cInverse(b));
}

vec2 cExp(vec2 z) {
	return vec2(exp(z.x) * cos(z.y), exp(z.x) * sin(z.y));
}

vec2 cLog(vec2 a) {
	float b =  atan(a.y,a.x);
	if (b>0.0) b-=2.0*3.1415;
	return vec2(log(length(a)),b);
}

vec2 cSqr(vec2 z) {
	return vec2(z.x*z.x-z.y*z.y,2.*z.x*z.y);
}

vec2 cSin(vec2 z) {
  return vec2(sin(z.x)*cosh(z.y), cos(z.x)*sinh(z.y));
}

vec2 cCos(vec2 z) {
  return vec2(cos(z.x)*cosh(z.y), -sin(z.x)*sinh(z.y));
}

#define CA vec2(par[2].xy) //uniform vec2 CA; //slider[(-1,-1),(0,0),(1,1)]
#define CB vec2(par[3].xy) //uniform vec2 CB; //slider[(-1,-1),(0,0),(1,1)]
#define CC vec2(par[4].xy) //uniform vec2 CC; //slider[(-1,-1),(0,0),(1,1)]
#define CD vec2(par[5].xy) //uniform vec2 CD; //slider[(-1,-1),(0,0),(1,1)]

// 'Ducks' fractal by Samuel Monnier
// (Implementation by Syntopia)
// See http://www.algorithmic-worlds.net/blog/blog.php?Post=20110227
vec2 formula(vec2 z) {
	z = cSqr(z);
	z = cDiv((cMul(CA,z)+CB),(cMul(z,CC)+CD));
	return z;
}

// Number of iterations
#define Iterations iters //slider[1,200,1000]

// Skip this number of iterations before starting coloring
#define PreIterations color_iters //slider[0,1,100]
#define EscapeSize par[0].x //slider[0,5,11]
float escape = pow(10.0,EscapeSize);

#define ColR par[14].x // {min=0 max=4 step=.01}
#define ColG par[14].y // {min=0 max=4 step=.01}
#define ColB par[14].z // {min=0 max=4 step=.01}
#define ColC par[15].x // {min=0 max=2 step=.01}

vec3 c2d(vec2 c) {
	vec2 z = c;
	int i = 0;
	float ci = 0.0;
	float mean = 0.0;
	for (i = 0; i < Iterations; i++) {
		z = formula(z);
		if (i>PreIterations) mean+=length(z);
		if (dot(z,z)> escape) break;
	}
	mean/=float(i-PreIterations);
	ci =  1.0 - log2(.5*log2(mean/ColC));
	return vec3( .5+.5*cos(6.*ci+ColR),.5+.5*cos(6.*ci + ColG),.5+.5*cos(6.*ci +ColB) );
}

// multiply by 2.0 to get z=-1 plane projection.
vec2 stereographicSphereToPlane(vec3 p) {
	float n = dot(p,p)+1.0;
	return 2.0*vec2(p.x/(1.0-p.z),p.y/(1.0-p.z));
}

// for z = 0 projection.
vec3 stereographicPlaneToSphere(vec2 p) {
	float n = dot(p,p)+1.0;
	return vec3(2.0*p.xy,n-2.0)/n;
}

vec3 c2dx(vec2 p) {
	if (mod(length(p.x),1.0)<0.02) return vec3(0.0);
	if (mod(length(p),1.0)<0.02) return vec3(0.0);
	return vec3(p,1.0);
}

float d_sphere(vec3 p) {
	float de = (length(p-vec3(0.0,0.0,0.0))-1.0); // A sphere
	return de;
}

float d_floor(vec3 p) {
	return abs(p.z + 1.0);
}

float  d(vec3 p) {
	return min(d_sphere(p), d_floor(p));
}

vec3 c(vec3 z) {
	if (d_floor(z) < d_sphere(z)) {
		return c2d(z.xy);
	} else {
		return c2d(stereographicSphereToPlane(z));
	}
}

const float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec2 Eps = vec2(0, max(d_pos, normal_eps));
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

vec3 background_color(in vec3 vp) {
  // Also draw duck on dome..
  return c2d(stereographicSphereToPlane(normalize(vp)))*.7;
}

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

  //float falloff = pow(1.0 + l1_dist, 3.0);
  //vec3 c = vec3(0.3,0.3,0.4) * clamp(2.0 / falloff, 0.0, 1.0);
  vec3 c = blinn_phong(n, -dp, dir, col, specularColor);
  if (totalD < l1_dist + fudge) c = c * .5;
  return c;
}

vec3 getLight1() {
  vec3 l = L1_Vector;
  return l * vec3(cos(time),sin(time), -1.5);  // circle around
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
    clamp(pow(1.5*(L1_Size-d)/L1_Size, 3.0), 0.0, 0.98));
}

// Get base color at p, plus Blinn_phing and ambient occulusion.
vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist, float side, float m_zoom) {
  vec3 col = c(p);
  col = blinn_phong(n, -dp, normalize(vec3(1.0,0.6,0.7)+dp), col, specularColor);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist, side));
  col = shade(p, getLight1(), side, m_zoom, n, col, dp);
  return col;
}

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 30x speed.
void setup_stereo(INOUT(vec3,eye_in), INOUT(vec3,dp)) {
#if !defined(ST_NONE)
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * abs(speed), 0, 0, 0));
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

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){ return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

void main() {
  vec3 eye_in, dp; setup_stereo(eye_in, dp);

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
    if (ds < d_floor(p)) {
      sphereHit = true;
      n = normalize(p);  // easy accurate normal for sphere hit
    } else {
      n = vec3(0,0,1);
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
          !sphereHit &&
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
      if (ds < d_floor(p)) {
        sphereHit = true;
        n = normalize(p);  // easy accurate normal for sphere hit
      } else {
        n = vec3(0,0,1);
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

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(firstD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(finalCol, depth);
//  gl_FragColor = vec4(vec3(pow(finalCol.x, 1.0/2.2),pow(finalCol.y, 1.0/2.2), pow(finalCol.z, 1.0/2.2)), depth);
}

