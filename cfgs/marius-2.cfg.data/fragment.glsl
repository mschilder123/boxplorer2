// Mandelbox shader by Rrrola
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal

#define P0 p0                    // standard Mandelbox
//#define P0 vec4(par[1].x,par[1].y,par[2].y,1)  // Mandelbox Julia

#define SCALE par[0].y
#define MINRAD2 par[0].x

#define DIST_MULTIPLIER 0.5
#define MAX_DIST 4.0

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
vec3 backgroundColor = //vec3(0.07, 0.06, 0.16),
  vec3(.2, .2, .2),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(.9, 0.8, 0.4),
  glowColor = vec3(0.03, 0.9, 0.04),
  aoColor = vec3(0, 0, 0);

// precomputed constants
float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
vec4 scale = vec4(SCALE, SCALE, SCALE, abs(SCALE)) / minRad2;
float absScalem1 = abs(SCALE - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(SCALE), float(1-iters));

// Rotate around vector
float Angle = par[3].x;
float csat = cos(Angle);
float ssat = sin(Angle);
float usat = 1.-cos(Angle);
vec3 p0 = normalize(vec3(par[2].x, par[2].y, par[3].y));
mat3 RotationMatrix = mat3( p0.x*p0.x*usat + csat,      p0.x*p0.y*usat + p0.z*ssat, p0.x*p0.z*usat - p0.y*ssat,
                            p0.y*p0.x*usat - p0.z*ssat, p0.y*p0.y*usat + csat,      p0.y*p0.z*usat + p0.x*ssat,
                            p0.z*p0.x*usat + p0.y*ssat, p0.z*p0.y*usat - p0.x*ssat, p0.z*p0.z*usat + csat
                      );


// Compute the distance from `pos` to the Mandelbox.
float de_box(vec3 pos) {
  vec4 p = vec4(pos,1), p0 = p;  // p.w is the distance estimate

  for (int i=0; i<iters; i++) {
    // box folding: if (p>1) p = 2-p; else if (p<-1) p = -2-p;
//    p.xyz = abs(1.0+p.xyz) - p.xyz - abs(1.0-p.xyz);  // add;add;abs.add;abs.add (130.4%)
//    p.xyz = clamp(p.xyz*0.5+0.5, 0.0, 1.0) * 4.0 - 2.0 - p.xyz;  // mad.sat;mad;add (102.3%)
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;  // min;max;mad

    // sphere folding: if (r2 < minRad2) p /= minRad2; else if (r2 < 1.0) p /= r2;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);  // dp3,div,max.sat,mul
  
    // scale, translate
    p = p*scale + P0;
  p.xyz *= RotationMatrix;
  }
  return ((length(p.xyz) - absScalem1) / p.w - AbsScaleRaisedTo1mIters) * DIST_MULTIPLIER;
}

// parameterized by iters
float de_box_iters(vec3 pos, int niters) {
  vec4 p = vec4(pos,1), p0 = p;  // p.w is the distance estimate

  for (int i=0; i<niters; i++) {
    // box folding: if (p>1) p = 2-p; else if (p<-1) p = -2-p;
//    p.xyz = abs(1.0+p.xyz) - p.xyz - abs(1.0-p.xyz);  // add;add;abs.add;abs.add (130.4%)
//    p.xyz = clamp(p.xyz*0.5+0.5, 0.0, 1.0) * 4.0 - 2.0 - p.xyz;  // mad.sat;mad;add (102.3%)
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;  // min;max;mad

    // sphere folding: if (r2 < minRad2) p /= minRad2; else if (r2 < 1.0) p /= r2;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);  // dp3,div,max.sat,mul
  
    // scale, translate
    p = p*scale + P0;
  p.xyz *= RotationMatrix;
  }
  return ((length(p.xyz) - absScalem1) / p.w - pow(abs(SCALE), float(1-niters))) * DIST_MULTIPLIER;
}

// Compute the color at `pos`.
vec3 color_box_iters(vec3 pos, int niters) {
  vec3 p = pos, p0 = p;
  float trap = 1.0;

  for (int i = 0; i < niters; i++) {
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale.xyz + P0.xyz;
//  p.xyz *= RotationMatrix;   // comment out for cool color interference, when rotated.
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}


//float normal_eps = 0.00001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos, float normal_eps, int niters) {
  vec4 Eps = vec4(0, normal_eps, 2.0*normal_eps, 3.0*normal_eps);
  return normalize(vec3(
  // 2-tap forward differences, error = O(eps)
//    -d_pos+de_box(pos+Eps.yxx),
//    -d_pos+de_box(pos+Eps.xyx),
//    -d_pos+de_box(pos+Eps.xxy)

  // 3-tap central differences, error = O(eps^2)
    -de_box_iters(pos-Eps.yxx, niters)+de_box_iters(pos+Eps.yxx, niters),
    -de_box_iters(pos-Eps.xyx, niters)+de_box_iters(pos+Eps.xyx, niters),
    -de_box_iters(pos-Eps.xxy, niters)+de_box_iters(pos+Eps.xxy, niters)

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
vec3 blinn_phong(vec3 normal, vec3 view, vec3 light, vec3 diffuseColor, float dist) {
  vec3 halfLV = normalize(light + view);
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return clamp(dif*dist,0.,1.)*diffuseColor + (spe/dist)*specularColor;
}


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n, int niters) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = de_box_iters(p + n*dist, iters);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_eps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

void main() {
  // Interlaced stereoscopic eye fiddling
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);

  vec3 p = eye_in, dp = normalize(dir);

  float totalD = 0.0, D = 0.0;

  int ncolor_iters = color_iters;

  // Intersect the view ray with the Mandelbox using raymarching.

  int steps;
  float m_dist = min_dist;

  for (steps = 0; steps < max_steps; steps++) {
    D = de_box(p + totalD * dp);
    if (D < m_dist) break;
    totalD += D;
  if (D > MAX_DIST) break;
  m_dist = min_dist * pow(totalD * par[8].x + 1.0, 100.*par[8].y);  // Vary LoD
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  float delta = m_dist / min_dist;
  int niters = int(clamp(float(iters) / pow(delta,par[9].x), 11., float(iters)));

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, D, m_dist * 10., iters);  // should look around in screen-space?

    col = color_box_iters(p, ncolor_iters);
#if 0
    int(clamp(
      float(ncolor_iters)/(2.*delta),
      5.0,
      float(color_iters))));
#endif

    col = blinn_phong(n, -dp, normalize(eye_in+vec3(0.,1.,0.)+dp), col, 1. /*clamp(delta/5., 1., 5.)*/);
    col = mix(aoColor, col, clamp(ambient_occlusion(p, n, niters), 0., 1.));

  // fog
  col = mix(col, backgroundColor, clamp(pow(1.+totalD*par[9].y, 2.)-1., 0., 1.));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor, clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
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
