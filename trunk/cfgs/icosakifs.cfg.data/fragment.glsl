// Mandelbox shader by Rrrola
// Original formula by Tglad
// - http://www.fractalforums.com/3d-fractal-generation/amazing-fractal
//bermarte:k-icosa,octa,dodeca formulas by knighty
#define d de_icosa

#define P0 p0                    // standard Mandelbox
//#define P0 vec4(par[1].x,par[1].y,par[2].y,1)  // Mandelbox Julia

#define SCALE par[0].y// {min=-3 max=3 step=.001}
#define MINRAD2 par[0].x// {min=0 max=3 step=.001}

#define CX par[6].x// {min=-3 max=3 step=.001}
#define CY par[7].x// {min=-3 max=3 step=.001}
//#define CZ par[8].x// {min=-3 max=3 step=.001}

#define prefolds par[5].x// {min=0 max=5 step=1}

#define Angle par[3].x// {min=-3 max=6 step=.001}
#define scale par[4].x// {min=0 max=9 step=.001}

#define _PHI_  (0.5*(1.+sqrt(5.)))//GOLDEN RATIO
#define ftris normalize(vec3(_PHI_,1.,0.))//VERTEX> ftris[3]
#define ftria normalize(vec3(_PHI_,0.,0.))//middle of edge> ftria[3]
#define ftric normalize(vec3( 1./3.*(1.+2.*_PHI_),0.,_PHI_/3.))//center icosa ftric[3]

#define IN3  (1./sqrt(14.+6.*sqrt(5.)))
#define n3 normalize(vec3(IN3*_PHI_,-IN3*(pow(_PHI_,2.)),-IN3*(2.*_PHI_+1.)))

#define DIST_MULTIPLIER 1.0
#define MAX_DIST 4.0

// Camera position and direction.
varying vec3 eye, dir;
varying float zoom;
uniform float xres, yres, time, speed;

// Interactive parameters.
uniform vec3 par[10];

uniform float
  min_dist,           // Distance at which raymarching stops.
  ao_eps,             // Base distance at which ambient occlusion is estimated.
  ao_strength,        // Strength of ambient occlusion.
  glow_strength,      // How much glow is applied after max_steps.
  dist_to_color;      // How is background mixed with the surface color after max_steps.

uniform int iters,    // Number of fractal iterations.
  color_iters,        // Number of fractal iterations for coloring.
  max_steps,          // Maximum raymarching steps.
  frameno;

// Colors. Can be negative or >1 for interestiong effects.
vec3 backgroundColor = vec3(.07, 0.06, 0.16),
  surfaceColor1 = vec3(0.95, 0.64, 0.1),
  surfaceColor2 = vec3(0.89, 0.95, 0.75),
  surfaceColor3 = vec3(0.55, 0.06, 0.03),
  specularColor = vec3(1.0, 0.8, 0.4),
  glowColor = vec3(0.03, 0.09, 0.4),
  aoColor = vec3(0, 0, 0);

// precomputed constants
float minRad2 = clamp(MINRAD2, 1.0e-9, 1.0);
vec4 scale2 = vec4(SCALE, SCALE, SCALE, abs(SCALE)) / minRad2;
float absScalem1 = abs(SCALE - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(SCALE), float(1-iters));

// rotations
//float Angle = par[3].x;
float csat = cos(Angle);
float ssat = sin(Angle);
float usat = 1.-cos(Angle);
vec3 z0 = normalize(vec3(par[2].x, par[2].y, par[3].y));
mat3 RotationMatrix = mat3(
z0.x*z0.x*usat + csat, z0.x*z0.y*usat + z0.z*ssat, z0.x*z0.z*usat - z0.y*ssat,
z0.y*z0.x*usat - z0.z*ssat, z0.y*z0.y*usat + csat, z0.y*z0.z*usat + z0.x*ssat,
z0.z*z0.x*usat + z0.y*ssat, z0.z*z0.y*usat - z0.x*ssat, z0.z*z0.z*usat + csat);

float de_icosa(vec3 z0) {
  int i; float t=0.0; float r=0.0;
    
  for(i=0;i<int(prefolds);i++) { //Pre-fold maxFoldIterations
    z0.yz = abs(z0.yz);

    t = dot(z0, n3);
    if(t<0.0){z0.x-=2.*t*n3[0];z0.y-=2.*t*n3[1];z0.z-=2.*t*n3[2];}
  }
  for (i=0;i<iters && r<20.0;i++) {
    z0 = z0*RotationMatrix;

    z0 = abs(z0);
    
    t = dot(z0, n3);
    if(t<0.){z0.x-=2.*t*n3[0];z0.y-=2.*t*n3[1];z0.z-=2.*t*n3[2];}

    z0.yz = abs(z0.yz);

    t = dot(z0, n3);
    if(t<0.){z0.x-=2.*t*n3[0];z0.y-=2.*t*n3[1];z0.z-=2.*t*n3[2];}

    z0.yz = abs(z0.yz);

    t = dot(z0, n3);
    if(t<0.){z0.x-=2.*t*n3[0];z0.y-=2.*t*n3[1];z0.z-=2.*t*n3[2];}

    // Stretching:
    // for an icosahedron: stc[]=ftris[];
    // for a dodecahedron: stc[]=ftric[];
    // In general you can choose: stc[]=a*ftris[]+b*ftria[]+c*ftric[]; where: a+b+c=1;
    z0.x=scale*z0.x-ftris[0]*(scale-CX);
    z0.y=scale*z0.y-ftris[1]*(scale-CY);
    z0.z=scale*z0.z-ftris[2]*(scale);

    r = dot(z0, z0);
  }

  float k= (sqrt(r)-2.0)*pow(scale,float(-i));
  return k;
}

// Compute the color at `pos`.
vec3 color(vec3 pos) {
  vec3 p = pos, p0 = p;
  float trap = 1.0;

  for (int i=0; i<color_iters; i++) {
    p.xyz = clamp(p.xyz, -1.0, 1.0) * 2.0 - p.xyz;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(minRad2/r2, minRad2), 0.0, 1.0);
    p = p*scale2.xyz + P0.xyz;
    trap = min(trap, r2);
  }
  // c.x: log final distance (fractional iteration count)
  // c.y: spherical orbit trap at (0,0,0)
  vec2 c = clamp(vec2( 0.33*log(dot(p,p))-1.0, sqrt(trap) ), 0.0, 1.0);

  return mix(mix(surfaceColor1, surfaceColor2, c.y), surfaceColor3, c.x);
}


float normal_eps = 0.0001;

// Compute the normal at `pos`.
// `d_pos` is the previously computed distance at `pos` (for forward differences).
vec3 normal(vec3 pos, float d_pos) {
  vec2 Eps = vec2(0, d_pos);
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
  float spe = pow(max( dot(normal, halfLV), 0.0 ), 32.0);
  float dif = dot(normal, light) * 0.5 + 0.75;
  return dif*diffuseColor + spe*specularColor;
}


// Ambient occlusion approximation.
float ambient_occlusion(vec3 p, vec3 n) {
  float ao = 1.0, w = ao_strength/ao_eps;
  float dist = 2.0 * ao_eps;

  for (int i=0; i<5; i++) {
    float D = d(p + n*dist);
    ao -= (dist-abs(D)) * w;
    w *= 0.5;
    dist = dist*2.0 - ao_eps;  // 2,3,5,9,17
  }
  return clamp(ao, 0.0, 1.0);
}

// ytalinflusa's noise [0..1>
float pnoise(vec2 pt){return mod(pt.x*(pt.x+0.15731)*0.7892+pt.y*(pt.y+0.13763)*0.8547,1.0); }

uniform float focus;  // {min=-10 max=30 step=.1} Focal plane devation from 30x speed.
void setup_stereo(inout vec3 eye_in, inout vec3 dp) {
#if !defined(ST_NONE)
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 4.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed, 0, 0, 0));
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

void main() {
  vec3 eye_in, dp; setup_stereo(eye_in, dp);
  float m_zoom = zoom * 0.5 / xres;
  vec3 p = eye_in;

  float noise = pnoise(gl_FragCoord.xy);
  float D = d(p);
  float side = 1.0; //sign(D); interior not good w/ this DE..
  float totalD = side * D * noise;
  float m_dist = m_zoom * totalD;

  // Intersect the view ray with the fractal using raymarching.
  int steps;
  for (steps=0; steps<max_steps; steps++) {
    D = side * d(p + totalD * dp);

    if (D < m_dist) break;
    totalD += D;
    if (totalD > MAX_DIST) break;

    m_dist = m_zoom * totalD;
  }

  p += totalD * dp;

  // Color the surface with Blinn-Phong shading, ambient occlusion and glow.
  vec3 col = backgroundColor;

  // We've got a hit or we're not sure.
  if (D < MAX_DIST) {
    vec3 n = normal(p, m_dist);
    col = color(p);
    col = blinn_phong(n, -dp, normalize(eye_in+vec3(0,1,0)+dp), col);
    col = mix(aoColor, col, ambient_occlusion(p, n));

    // We've gone through all steps, but we haven't hit anything.
    // Mix in the background color.
    if (D > m_dist) {
      col = mix(col, backgroundColor,
                clamp(log(D/m_dist) * dist_to_color, 0.0, 1.0));
    }
  }

  // Glow is based on the number of steps.
  col = mix(col, glowColor, float(steps)/float(max_steps) * glow_strength);

  float zNear = abs(speed);
  float zFar = 65535.0 * zNear;
  float a = zFar / (zFar - zNear);
  float b = zFar * zNear / (zNear - zFar);
  float depth = (a + b / clamp(totalD/length(dir), zNear, zFar));
  gl_FragDepth = depth;
  gl_FragColor = vec4(col, depth);
}
