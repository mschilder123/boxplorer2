// Menger-sphere shader by marius. Parts by knighty, Rrrola.
// See http://www.fractalforums.com/3d-fractal-generation/revenge-of-the-half-eaten-menger-sponge/

#define MAX_DIST 4.0

// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;  // Distance at which raymarching stops. {min=1e-08 max=1e-04 step=1e-08}
uniform float ao_eps;  // Base distance at which ambient occlusion is estimated. {min=0 max=.001 step=.000001}
uniform float ao_strength;  // Strength of ambient occlusion. {min=0 max=.01 step=.0001}
uniform float glow_strength;  // How much glow is applied after max_steps. {min=0 max=10 step=.05}
uniform float dist_to_color;  // How is background mixed with the surface color after max_steps. {min=0 max=10 step=.05}

uniform float speed;  // {min=1e-06 max=.1 step=1e-06}

uniform int iters;  // Number of fractal iterations. {min=1 max=100}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=1 max=100}
uniform int max_steps;  // Maximum raymarching steps. {min=1 max=200}

// Colors. Can be negative or >1 for interesting effects.
#define backgroundColor par[7]
#define surfaceColor1 par[2]
//#define surfaceColor2 par[3]
//#define surfaceColor3 par[4]
#define specularColor par[5]
#define glowColor par[6]
#define lightVector par[1]

#define MOD par[8].x  // {min=1 max=5 step=.1}
#define SCALE par[8].y  // {min=1 max=5 step=.1}
#define INTRA_D par[8].z  // {min=.01 max=10.0 step=.01}
#define INITIAL_K par[0].x  // {min=0 max=5 step=.1}
#define DE_EPS par[0].y  // {min=.0000001 max=.001 step=.0000001}

vec3 aoColor = vec3(0, 0, 0);

//Fuctions to call.
#define d d_SphereSponge
#define color color_SphereSponge

float d_SphereSponge(vec3 pos) {
  float k = INITIAL_K;
  float d = -100.0;
  for(int i=0; i<iters; i++) {
    vec3 x = mod(pos * k, MOD) - .5 * MOD;
    float r = length(x);
    float d1 = (INTRA_D - r) / k;
    d = max(d, d1);
    k *= SCALE;
  }
  return d;
}

vec3 color_SphereSponge(vec3 pos) {
  return surfaceColor1;
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


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, float DistAtp, float side) {
  float ao_ed=DistAtp*100.0;
  float ao = 1.0, w = ao_strength/ao_ed;//ps;
  float dist = 2.0 * ao_ed;//ps;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_ed;//ps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);

  vec3 p = eye_in, dp = normalize(dir);

  float totalD = 0.0, D = 1000.0, extraD = 0.0, lastD;
  D = d(p + totalD * dp);
  float side = sign(D);
  D = D*side; D *= 0.95;
  totalD = D;

  // Intersect the view ray with the Mandelbox using raymarching.
  int steps;
  for (steps=0; steps<max_steps && D > max(DE_EPS*totalD,0.0000152587890625) && totalD < MAX_DIST; steps++) {
    D = side * d(p + totalD * dp);
    D = D - min_dist * (totalD + D) ;
    totalD += D;
  }

  p += totalD * dp;
  //normal_eps = max(normal_eps,1.0*D);
  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    D = max(D, max(DE_EPS*totalD,0.0000152587890625));
    vec3 n = side * normal(p, D);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(eye_in+lightVector+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, D, side));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > min_dist) {
      col = mix(col, backgroundColor, clamp(log(D/min_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  // Write zBuffer and pixel
  float zFar = 5.0;
  float zNear = 0.0001;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
