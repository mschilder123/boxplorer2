// Mandalabeth 1522 shader by marius.
// Original shader by Rrrola for mandelbox, part of boxplorer.

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
uniform float xres;
varying float zoom;

uniform int iters;  // Number of fractal iterations. {min=1 max=100}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=1 max=100}
uniform int max_steps;  // Maximum raymarching steps. {min=1 max=200}

#define SIGN par[0].z  // {min=-1 max=1 step=2}

// Colors. Can be negative or >1 for interesting effects.
#define backgroundColor par[7]
#define surfaceColor1 par[2]
#define surfaceColor2 par[3]
#define surfaceColor3 par[4]
#define specularColor par[5]
#define glowColor par[6]
#define lightVector par[1]

vec3 aoColor = vec3(0, 0, 0);

//Fuctions to call.
#define d d_MandalaBeth1522  //d_MandalaBeth1032
#define color color_MandalaBeth

float d_MandalaBeth1522(vec3 pos) {
   float u = pos.x, v = pos.y, w = pos.z;
   float sign = clamp(SIGN, -1., 1.);
   float t = sqrt(5.),
   x = u, x2 = x*x,
   y = v, y2 = y*y,
   z = w, z2 = z*z,
   r2 = x2+y2+z2,dr = 1.;
   for(int i = 0; (i < iters) && (r2 <= 4.); i++)
   {
      dr = 1.+2.5*r2*r2*dr;
      float
      x4 = x2*x2,
      y4 = y2*y2,
      z4 = z2*z2,
      uu = u-sign*5./32.*x*((5.+7.*t)*y4+(5.-7.*t)*z4-60.*y2*z2+2.*x2*((5.-7.*t)*y2+(5.+7.*t)*z2)-2.*x4),
      vv = v-sign*5./32.*y*((5.+7.*t)*z4+(5.-7.*t)*x4-60.*z2*x2+2.*y2*((5.-7.*t)*z2+(5.+7.*t)*x2)-2.*y4),
      ww = w-sign*5./32.*z*((5.+7.*t)*x4+(5.-7.*t)*y4-60.*x2*y2+2.*z2*((5.-7.*t)*x2+(5.+7.*t)*y2)-2.*z4);

      x = uu; x2 = x*x;
      y = vv; y2 = y*y;
      z = ww; z2 = z*z;
      r2 = x2+y2+z2;
   }

   float r = sqrt(r2);
   return .4*log(r)*r/(2.*dr);
}

float d_MandalaBeth1032(vec3 pos) {
  float u = pos.x, v = pos.y, w = pos.z;
  float sign = clamp(SIGN, -1., 1.);
  float t = sqrt(5.),
  x = u, x2 = x*x,
  y = v, y2 = y*y,
  z = w, z2 = z*z,
  r2 = x2+y2+z2,dr = 1.;

  for(int i = 0; (i < iters) && (r2 <= 4.); i++)
  {
     dr = 1.0+3.5*r2*r2*r2*dr;
     float
     x4 = x2*x2, x6 = x4*x2,
     y4 = y2*y2, y6 = y4*y2,
     z4 = z2*z2, z6 = z4*z2,
     uu = u+sign*56./81.*x*(2.*(4.+3.*t)*y6+2.*(4.-3.*t)*z6-15.*(5.-3.*t)*y4*z2-15.*(5.+3.*t)*y2*z4
                      -15.*x2*((1.+3.*t)*y4+(1.-3.*t)*z4-20.*y2*z2)
                      -3.*x4*((7.-9.*t)*y2+(7.+9.*t)*z2)
                      +2.*x6),
     vv = v+sign*56./81.*y*(2.*(4.+3.*t)*z6+2.*(4.-3.*t)*x6-15.*(5.-3.*t)*z4*x2-15.*(5.+3.*t)*z2*x4
                      -15.*y2*((1.+3.*t)*z4+(1.-3.*t)*x4-20.*z2*x2)
                      -3.*y4*((7.-9.*t)*z2+(7.+9.*t)*x2)
                      +2.*y6),
     ww = w+sign*56./81.*z*(2.*(4.+3.*t)*x6+2.*(4.-3.*t)*y6-15.*(5.-3.*t)*x4*y2-15.*(5.+3.*t)*x2*y4
                      -15.*z2*((1.+3.*t)*x4+(1.-3.*t)*y4-20.*x2*y2)
                      -3.*z4*((7.-9.*t)*x2+(7.+9.*t)*y2)
                      +2.*z6);
     x = uu; x2 = x*x;
     y = vv; y2 = y*y;
     z = ww; z2 = z*z;
     r2 = x2+y2+z2;
  }
  float r = sqrt(r2);
  return log(r)*r/(50.*dr);
}


// Compute the color at `pos`.
vec3 color_MandalaBeth(vec3 pos) {
   float u = pos.x, v = pos.y, w = pos.z;
   float sign = clamp(SIGN, -1., 1.);
   float t = sqrt(5.),
   x = u, x2 = x*x,
   y = v, y2 = y*y,
   z = w, z2 = z*z,
   r2;
   float trap = 1.0;   
   for(int i = 0; i < color_iters; i++)
   {
      float
      x4 = x2*x2,
      y4 = y2*y2,
      z4 = z2*z2,
      uu = u-sign*5./32.*x*((5.+7.*t)*y4+(5.-7.*t)*z4-60.*y2*z2+2.*x2*((5.-7.*t)*y2+(5.+7.*t)*z2)-2.*x4),
      vv = v-sign*5./32.*y*((5.+7.*t)*z4+(5.-7.*t)*x4-60.*z2*x2+2.*y2*((5.-7.*t)*z2+(5.+7.*t)*x2)-2.*y4),
      ww = w-sign*5./32.*z*((5.+7.*t)*x4+(5.-7.*t)*y4-60.*x2*y2+2.*z2*((5.-7.*t)*x2+(5.+7.*t)*y2)-2.*z4);

      x = uu; x2 = x*x;
      y = vv; y2 = y*y;
      z = ww; z2 = z*z;
      r2 = x2+y2+z2;
      trap = min(trap, r2);
   }

   //vec2 c = clamp(vec2( 5.*log(r2)-1.0, sqrt(trap)), 0., 6.);
   vec2 c = normalize(vec2(log(r2), sqrt(trap)));
   return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
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
float ambient_occlusion(vec3 p, vec3 n, float side) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = side * d(p + n*dist);
    ao -= (dist-D) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_eps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}


void main() {
  vec3 eye_in = eye;
  eye_in += 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed *
      vec3(gl_ModelViewMatrix[0]);

  vec3 p = eye_in, dp = normalize(dir);

  float totalD = 0.0, D = 3.4e38, extraD = 0.0, lastD;
  float cutD;

  vec3 col = backgroundColor;
  int steps = 0;
  float m_dist = min_dist;
  float m_zoom = .5 * zoom / xres;

  cutD = (p.z - par[0].x) / -dp.z;
  if (p.z < par[0].x) cutD = -1.;
  if (cutD > 0.) totalD = cutD;

  // Intersect the view ray with the fratal using raymarching.
  for (; steps<max_steps; steps++) {
    lastD = D;
    D = d(p + totalD * dp);
    if (extraD > 0.0 && D < extraD) {
        totalD -= extraD;
        extraD = 0.0;
        D = 3.4e38;
        steps--;
        continue;
    }
    if (D < m_dist) break;
    if (D > MAX_DIST) break;

    totalD += D;
    totalD += extraD = 0.096 * D * (D+extraD)/lastD;
	m_dist = m_zoom * totalD;
  }

float z = (p+totalD*dp).z;
if (z - par[0].x < -m_dist) {
  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.


  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, D);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(eye_in+vec3(lightVector.x,lightVector.y,-lightVector.z)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n, .5));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor, clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);
 } else {
  if (cutD>0. && D<MAX_DIST) {
    p+=cutD*dp;
    col = color(p.xyz);
    vec3 n = vec3(0,0,1.);  // fixed normal of cutting plane, no ambient occlusion..
    col = blinn_phong(n, -dp, normalize(eye_in+vec3(lightVector.x,lightVector.y,-lightVector.z)+dp), col);
  } else totalD = MAX_DIST;
 }

  // Write Z-buffer
  float zNear = abs(speed);
  float zFar = 65535.0*zNear;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
