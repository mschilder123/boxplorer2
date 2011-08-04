// menger and then some shader.
// Original shader by rrrola for mandelbox
// bermarte: formula from Knighty
// marius: refactored w/ reflections, background dome, ssponge, combi etc.
// marius: massaged so can be compiled as plain C++.

#define d de_menger // combi,menger,mandelbox,ssponge  // distance estimator
#define c c_menger  // color at position p

#define MAX_DIST 6.0
#define ULP 0.000000059604644775390625

#ifndef PI
#define PI 3.14159265
#endif
#ifndef INOUT
#define INOUT(a,b) inout a b
#endif

// Camera position and direction.
varying vec3 eye, dir;
varying float zoom;
uniform float xres;

// Interactive parameters.
uniform vec3 par[20];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform float speed;  // eye separation really.

uniform int iters,    // Number of fractal iterations.
  color_iters,        // Number of reflection rays.
  max_steps;          // Maximum raymarching steps.

// Colors. Can be negative or >1 for interesting effects.
#define backgroundColor par[10]
#define surfaceColor1 par[11]
#define surfaceColor2 par[12]
#define surfaceColor3 par[13]
#define specularColor par[14]
#define glowColor par[15]
#define aoColor par[16]

#define SHINE par[9].x  // {min=0 max=1 step=.01} reflection contribution drop-off
#define SIGMA par[9].y  //{min=0 max=1 step=.01} surface roughness
#define ALBEDO par[9].z  //{min=0 max=1 step=.01} surface reflection

float de_menger(vec3 z0) {
#define offsetVector par[2]
#define M_ITERS par[1].x //{min=0 max=20 step=1} iterations
#define M_SIZE par[1].y //{min=0 max=5 step=.01} size of menger box
  float scale=2.;  // size of holes
  float r=1.0;
  z0 /= M_SIZE;
  int n_iters = int(M_ITERS);
  for (int i=0;i<n_iters;i++){
    z0 = abs(z0) + par[2];
    if( z0.x < z0.y){z0 = z0.yxz;}
    if( z0.x < z0.z){z0 = z0.zyx;}
    if( z0.y < z0.z){z0 = z0.xzy;}
    r = z0.x;
    z0 += z0*scale - scale;
    if(z0.z<-0.5*scale) z0.z+=scale;
  }
  return (r-1.0)*pow(scale+1.f,float(1-n_iters))*M_SIZE;
}

vec3 c_menger(vec3 p) {
  return surfaceColor1;  // Boring but reflections make it interesting.
}

#define SCALE par[0].y  //{min=-3 max=3 step=.01}
#define MINRAD2 par[0].x  //{min=0 max=1 step=.001}


float de_mandelbox(vec3 pos) {
  float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
  vec4 scale = vec4(SCALE, SCALE, SCALE, abs(SCALE)) / minRad2;
  vec4 p = vec4(pos,1.0), p0 = p;  // p.w is the distance estimate
  for (int i=0; i<iters; i++) {
    p = vec4(clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz, p.w);
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
	p = p*scale + p0;
  }
  return ((length(p.xyz) - abs(SCALE - 1.0)) / p.w
           - pow(abs(SCALE), float(1-iters))) * .85;
}

// Infinite construction menger sphere sponge.
float de_ssponge(vec3 pos) {
#define MOD par[8].x  // {min=1 max=5 step=.01}
#define SSCALE par[8].y  // {min=1 max=5 step=.01}
#define INTRA_D par[8].z  // {min=.01 max=10.0 step=.01}
#define INITIAL_K par[0].z  // {min=0 max=5 step=.01}
  float k = INITIAL_K;
  float d = -100.0;
  for(int i=0; i<iters; i++) {
    vec3 x = mod(pos * k, MOD) - .5 * MOD;
    float r = length(x);
    float d1 = (INTRA_D - r) / k;
    d = max(d, d1);
    k *= SSCALE;
  }
  return d;
}

float de_combi(vec3 z0) {
  // Use min() for union, max() for intersection of shapes.
  return min(max(de_ssponge(z0),de_mandelbox(z0)), de_menger(z0));
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
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0f);
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
float ambient_occlusion(vec3 p, vec3 n, float totalD, float m_dist) {
  float ao_ed = totalD*ao_eps/m_dist;
  float ao = 1.0, w = ao_strength/ao_ed;
  float dist = 2.0 * ao_ed;

  for (int i=0; i<5; i++) {
    float D = d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

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

// Get base color at p, plus Blinn_phing and ambient occulusion.
vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist) {
  vec3 col = c(p);
  col = blinn_phong(n, -dp, normalize(vec3(0,1,0)+dp), col);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist));
  return col;
}

#ifndef _FAKE_GLSL_

// Intersect direction ray w/ large encapsulating sphere.
// Sphere map texture onto it.
uniform sampler2D bg_texture;
uniform int use_bg_texture;
vec3 background_color(in vec3 vp) {
#define BG_BLUR par[1].z  //{min=0 max=8 step=.1}
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
	return texture2DLod(bg_texture, vec2(u,v), BG_BLUR).xyz;
}

void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);

  vec3 p = eye_in, dp = normalize(dir);
  float m_zoom = length(dir) * zoom * .5 / xres;  // screen error at dist 1.

  float totalD = d(p);
  float side = sign(totalD);
  totalD *= side * .5;  // start with conservative step.
  float m_dist = max(min_dist, m_zoom * totalD);

  int steps = 0;  // number of marching steps we've taken for this ray.
  float colFactor = SHINE;

  // March first ray.
  steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);
  vec3 rayCol, n;
  if (totalD < MAX_DIST) {
    p += totalD * dp;
	n = normal(p, m_dist * .5);
	rayCol = rayColor(p, dp, n, totalD, m_dist);
	rayCol = mix(rayCol, glowColor, float(steps)/float(max_steps) * glow_strength);
  } else {
    rayCol = background_color(dp);
  }

  float firstD = totalD;
  vec3 finalCol = rayCol;

  // March reflected ray a couple of times.
  for (int ray = 1; ray < color_iters && totalD < MAX_DIST; ++ray) {
    dp = reflect(dp, n);  // reflect view direction
    p += (-totalD + 2.0 * m_dist) * dp;  // reproject eye

    steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);
	if (totalD < MAX_DIST) {
      p += totalD * dp;
      n = normal(p, m_dist * .5);
	  rayCol = rayColor(p, dp, n, totalD, m_dist);
	  rayCol = mix(rayCol, glowColor, float(steps)/float(max_steps) * glow_strength);
    } else {
	  rayCol = background_color(dp);
	}

	finalCol = mix_reflection(n, -dp, finalCol, rayCol, colFactor);
	colFactor *= SHINE * SHINE;  // reflection drop-off.
  }

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(firstD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(finalCol, depth);
}

#endif  // _FAKE_GLSL_