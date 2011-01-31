// Mandelbox shader by Rrrola
// Original formula by Tglad
//mandelbulb mod.by visual/bermarte>thanx 2 Xyrus and 2 subblue
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal

#define DE de_bulb
#define COLOR color_bulb

#define Power par[0].x  // {min=0 max=15 step=.5}
#define Bailout par[0].y  // {min=0 max=100 step=1}

#define MAX_DIST 4.0

#define LoDbase par[8].x  // {min=1 max=20 step=.25}
#define LoDpow par[8].y  // {min=0 max=50 step=.001}

#define JuliaFactor par[0].z
#define JuliaVector par[2]

// Camera position and direction.
varying vec3 eye, dir;

// Interactive parameters.
uniform vec3 par[10];

uniform float min_dist;  // Distance at which raymarching stops. {min=1e-08 max=1e-03 step=1e-08}
uniform float ao_eps;  // Base distance at which ambient occlusion is estimated. {min=0 max=.001 step=.000001}
uniform float ao_strength;  // Strength of ambient occlusion. {min=0 max=.01 step=.0001}
uniform float glow_strength;  // How much glow is applied after max_steps. {min=0 max=10 step=.05}
uniform float dist_to_color;  // How is background mixed with the surface color after max_steps. {min=0 max=10 step=.05}

uniform float speed;  // {min=1e-06 max=.1 step=1e-06}

uniform int iters;  // Number of fractal iterations. {min=1 max=100}
uniform int color_iters;  // Number of fractal iterations for coloring. {min=1 max=100}
uniform int max_steps;  // Maximum raymarching steps. {min=1 max=200}

// Colors. Can be negative or >1 for interestiong effects.
#define backgroundColor par[3]
#define surfaceColor1 par[4]
#define surfaceColor2 par[5]
#define surfaceColor3 par[6]
#define specularColor par[7]
#define glowColor par[1]

#define lightVector par[9]

vec3  aoColor = vec3(.1, .1, .1);

// Compute the distance from `pos` to the Mandelbulb.
float de_bulb(vec3 z0) {
	vec3 c = mix(z0, JuliaVector, JuliaFactor);
	vec3 z=z0;

    	float pd = Power - 1.0;		

	float r	 = length(z);
	float th = atan(z.x, z.y);
	float ph = asin(z.z / r);

	vec3 dz = vec3(0,0,0);
	float ph_dz = 0.0;
	float th_dz = 0.0;
	float r_dz = 1.0;
	float powR, powRsin;

	for (int n = 0; n < iters; n++) {
		powR = Power * pow(r, pd);
		powRsin = powR * r_dz * sin(ph_dz + pd*ph);
		dz.x = powRsin * cos(th_dz + pd*th) + 1.0;
		dz.y = powRsin * sin(th_dz + pd*th);
		dz.z = powR * r_dz * cos(ph_dz + pd*ph);

		r_dz  = length(dz);
		th_dz = atan(dz.x, dz.y);
		ph_dz = acos(dz.z / r_dz);
		powR = pow(r, Power);
		powRsin = sin(Power*ph);
		z.x = powR * powRsin * cos(Power*th);
		z.y = powR * powRsin * sin(Power*th);
		z.z = powR * cos(Power*ph);
		z += c;

		r  = length(z);
		if (r > Bailout) break;

		th = atan(z.x, z.y)+10.0;
		ph = acos(z.z / r)+10.0;
	}

	float bulb=(0.5 * r * log(r)/r_dz);
	return bulb;
}

vec3 color_bulb(vec3 z0) {
	vec3 c = mix(z0, JuliaVector, JuliaFactor);
	vec3 z=z0;

    	float pd = Power - 1.0;		

	float r	 = length(z);
	float th = atan(z.x, z.y);
	float ph = asin(z.z / r);

	vec3 dz = vec3(0,0,0);
	float ph_dz = 0.0;
	float th_dz = 0.0;
	float r_dz = 1.0;
	float powR, powRsin;
	float trap = 1.0;

	for (int n = 0; n < color_iters; n++) {
		powR = Power * pow(r, pd);
		powRsin = powR * r_dz * sin(ph_dz + pd*ph);
		dz.x = powRsin * cos(th_dz + pd*th) + 1.0;
		dz.y = powRsin * sin(th_dz + pd*th);
		dz.z = powR * r_dz * cos(ph_dz + pd*ph);

		r_dz  = length(dz);
		th_dz = atan(dz.x, dz.y);
		ph_dz = acos(dz.z / r_dz);
		powR = pow(r, Power);
		powRsin = sin(Power*ph);
		z.x = powR * powRsin * cos(Power*th);
		z.y = powR * powRsin * sin(Power*th);
		z.z = powR * cos(Power*ph);
		z += c;

		r  = length(z);

		th = atan(z.x, z.y)+10.0;
		ph = acos(z.z / r)+10.0;
		trap = min(trap, r);
	}

	vec2 cl = clamp(vec2( .33*log(r)-1.0, trap), 0.0, 1.0);
	return mix(mix(surfaceColor1, surfaceColor2, cl.y), surfaceColor3, cl.x);
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
    -DE(pos-Eps.yxx)+DE(pos+Eps.yxx),
    -DE(pos-Eps.xyx)+DE(pos+Eps.xyx),
    -DE(pos-Eps.xxy)+DE(pos+Eps.xxy)

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
float ambient_occlusion(vec3 p, vec3 n) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = DE(p + n*dist);
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

  float odd = fract(gl_FragCoord.y * .5);
  float displace = (4. * odd - 1.) * speed;
  p += displace * vec3(gl_ModelViewMatrix[0]);

  float totalD = 0.0, D = 0.0;

  // Intersect the view ray with the Mandelbox using raymarching.

  int steps;
  float m_dist = min_dist;

  for (steps = 0; steps < max_steps; steps++) {
    D = DE(p + totalD * dp);
    if (D < m_dist) break;
    totalD += D;
    if (D > MAX_DIST) break;
    m_dist = min_dist * pow(totalD * LoDbase + 1.0, LoDpow);  // Vary LoD
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, D);
    col = COLOR(p);
    col = blinn_phong(n, -dp, normalize(eye_in+lightVector+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n));

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
