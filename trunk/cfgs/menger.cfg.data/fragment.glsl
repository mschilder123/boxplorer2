// menger shader.
// Original shader by rrrola for mandelbox
// bermarte: formula from Knigthy
// marius: refactored w/ reflections

#define d de_menger

#define MAX_DIST 4.0
#define ULP 0.000000059604644775390625
#define PI 3.14159265

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

#define SHINE par[9].x  // {min=0 max=1 step=.01}

float de_menger(vec3 z0) {
  //menger box
  //from Master-Knighty
  //brutal test
  int i;
  float scale=3.;
  float r=0.;
  for (i=0;i<iters;i++){
	vec3 zz0;
	z0.x=abs(z0.x);z0.y=abs(z0.y);z0.z=abs(z0.z);
	if( z0.x- z0.y<0.0){zz0.x=z0.y;z0.y=z0.x;z0.x=zz0.x;}
	if( z0.x- z0.z<0.0){zz0.x=z0.z;z0.z=z0.x;z0.x=zz0.x;}
	if( z0.y- z0.z<0.0){zz0.y=z0.z;z0.z=z0.y;z0.y=zz0.y;}
	zz0.x=z0.x-1.0;zz0.y=z0.y-1.0;zz0.z=z0.z-1.0;
	r=max(zz0.x,max(zz0.y,zz0.z));
	z0.x=z0.x*scale-1.0*(scale-1.0);
	z0.y=z0.y*scale-1.0*(scale-1.0);
	z0.z=z0.z*scale;
	if(z0.z>0.5*(scale-1.0)) z0.z-=(scale-1.0);
	}
  return r*pow(scale,1.0-float(i));
}

float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  //vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  vec2 Eps = vec2(0, max(d_pos, normal_eps));
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
  //  -2.0*d(pos-Eps.yxx)-3.0*d_pos+6.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
  //  -2.0*d(pos-Eps.xyx)-3.0*d_pos+6.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
  //  -2.0*d(pos-Eps.xxy)-3.0*d_pos+6.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)

  // 5-tap central differences, error = O(eps^4)
//    d(pos-Eps.zxx)-8.0*d(pos-Eps.yxx)+8.0*d(pos+Eps.yxx)-d(pos+Eps.zxx),
//    d(pos-Eps.xzx)-8.0*d(pos-Eps.xyx)+8.0*d(pos+Eps.xyx)-d(pos+Eps.xzx),
//    d(pos-Eps.xxz)-8.0*d(pos-Eps.xxy)+8.0*d(pos+Eps.xxy)-d(pos+Eps.xxz)
  ));
}

// Blinn-Phong shading model with rim lighting (diffuse light bleeding to the other side).
// `normal`, `view` and `light` should be normalized.
vec3 blinn_phong(in vec3 normal, in vec3 view, in vec3 light, in vec3 diffuseColor) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
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
int rayMarch(in vec3 p, in vec3 dp, inout float totalD, in float side, inout float m_dist, in float m_zoom) {
  int steps;
  for (steps = 0; steps < max_steps; ++steps) {
    float D = (side * d(p + totalD * dp) - totalD * m_dist) / (1.0 + m_dist);
    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;
    m_dist = max(2.*ULP, m_zoom * totalD);
  }
  return steps;
}

vec3 rayColor(vec3 p, vec3 dp, vec3 n, float totalD, float m_dist) {
  vec3 col = surfaceColor1;
  col = blinn_phong(n, -dp, normalize(vec3(0,1,0)+dp), col);
  col = mix(aoColor, col, ambient_occlusion(p, n, totalD, m_dist));
  return col;
}

// Intersect ray w/ large encapsulating sphere.
// If hit (always does), sphere map texture onto it.
uniform sampler2D bg_texture;
uniform int use_bg_texture;
vec3 background_color(in vec3 eye, in vec3 dir) {
  if (use_bg_texture == 0) return backgroundColor;
  float v = dot(eye, dir);
  float r = 100.0; // universe radius
  float d = r*r - (dot(eye, eye) - v*v);
  if (d < 0.0) {
    return backgroundColor;
  } else {
    // Got intersection, compute texture coords for that location.
	vec3 vp = normalize(eye + (v - sqrt(d)) * dir);
	vec3 vn = vec3(0.0, 1.0, 0.0);
	vec3 ve = vec3(1.0, 0.0, 0.0);
	float phi = acos(-dot(vn, vp));
	float v = phi / PI;
	float theta = (acos(dot(vp, ve) / sin(phi))) / (2. * PI);
	float u;
	if (dot(cross(vn, ve), vp) > 0.) {
	  u = theta;
	} else {
	  u = 1. - theta;
	}
	return texture2DLod(bg_texture, vec2(u,v),
	                    0.  // pick one mipmap down to kill seam?
                        ).xyz;
  }
}

void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);

  vec3 p = eye_in, dp = normalize(dir);
  float m_zoom = length(dir) * zoom * .5 / xres;  // screen error at dist 1.

  float D = d(p);
  float side = sign(D);
  float totalD = abs(D) * .5;
  float m_dist = m_zoom * totalD;

  int steps = 0;  // number of marching steps we've taken for this ray.
  vec3 finalCol = vec3(0.,0.,0.);
  float colFactor = 1.;
  float firstD = MAX_DIST;

  for (int ray = 0; ray < color_iters && totalD < MAX_DIST; ++ray) {
    steps += rayMarch(p, dp, totalD, side, m_dist, m_zoom);
	vec3 rayCol;
	if (totalD < MAX_DIST) {
      p += totalD * dp;
      vec3 n = normal(p, m_dist * .5);
	  rayCol = rayColor(p, dp, n, totalD, m_dist);
	  rayCol = mix(rayCol, glowColor, float(steps)/float(max_steps) * glow_strength);
      
      dp = reflect(dp, n);  // reflect direction
	  p += (-totalD + 2.0 * m_dist) * dp;  // reproject eye
    } else {
	  rayCol = background_color(p, dp);
	}
	finalCol = mix(rayCol, finalCol, 1. - colFactor);
	colFactor *= SHINE;
	firstD = min(totalD, firstD);  // track first dist for z-buffer
  }

  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(firstD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(finalCol, depth);
}
