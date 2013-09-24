//from http://glsl.heroku.com/e#9837.0 by @301z
//Hacked up by marius to work /w boxplorer2

uniform float xres, yres, speed, time;
varying vec3 eye, dir;  

uniform sampler2D bg_texture;
uniform int bg_weight;

uniform vec3 par[20];

#include "setup.inc"

#line 14
vec2 iResolution = vec2(xres, yres);

// by @301z

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Sphere {
  vec4 positionRadius;
  vec4 colorEmission;
  int reflectionType;
};

struct Hit {
  Sphere sphere;
  float distance;
};

const float kTwoPI = 2.0 * 3.14159;
const float kInf = 1e20;

const int kDiffuse = 0;
const int kSpecular = 1;
const int kRefraction = 3;

const Sphere kLeft = Sphere(vec4(1e4 + 1.0, 40.8, 81.6, 1e4), vec4(0.75, 0.25, 0.25, 0.0), kSpecular);
const Sphere kRight = Sphere(vec4(-1e4 + 99.0, 40.8, 81.6, 1e4), vec4(0.25, 0.25, 0.75, 0.0), kDiffuse);
const Sphere kBack = Sphere(vec4(50.0, 40.8, 1e4, 1e4), vec4(0.75, 0.75, 0.75, 0.0), kDiffuse);
const Sphere kFront = Sphere(vec4(50.0, 40.8, -1e4 + 170.0, 1e4), vec4(0.25, 0.75, 0.0, 0.0), kSpecular);
const Sphere kBottom = Sphere(vec4(50.0, 1e4, 81.6, 1e4), vec4(0.75, 0.75, 0.75, 0.0), kDiffuse);
const Sphere kTop = Sphere(vec4(50.0, -1e4 + 81.6, 81.6, 1e4), vec4(0.75, 0.75, 0.75, 0.0), kDiffuse);
const Sphere kMirror = Sphere(vec4(27.0, 16.5, 47.0, 16.5), vec4(0.999, 0.999, 0.999, 0.0), kSpecular);
const Sphere kGlass = Sphere(vec4(73.0, 16.5, 78.0, 16.5), vec4(0.999, 0.999, 0.999, 0.0), kRefraction);
// Note this is the only object that emits light (12.0)
const Sphere kLight = Sphere(vec4(50.0, 81.6 - 15.0, 81.6, 7.0), vec4(0.0, 0.0, 0.0, 12.0), kDiffuse);

Hit intersect(Hit hit, Ray ray, Sphere sphere) {
  vec3 op = sphere.positionRadius.xyz - ray.origin;
  float b = dot(op, ray.direction);
  float d = b * b - dot(op, op) + sphere.positionRadius.w * sphere.positionRadius.w;
  float ds = sqrt(d);
  Hit newHit = Hit(sphere, b + sign(ds - b) * ds);
  return all(bvec3(d > 0.0, newHit.distance > 0.01, newHit.distance < hit.distance)) ? newHit : hit;
}

Hit intersect(Ray ray) {
  Hit hit = Hit(kLight, kInf);
  hit = intersect(intersect(intersect(hit, ray, kLeft), ray, kRight), ray, kBack);
  hit = intersect(intersect(intersect(hit, ray, kFront), ray, kBottom), ray, kTop);
  hit = intersect(intersect(intersect(hit, ray, kMirror), ray, kGlass), ray, kLight);
  return hit;
}

vec3 cosineWeightedDirection(float u, float v, vec3 normal) {
  float r = sqrt(u);
  float angle = kTwoPI * v;
  vec3 sdir = cross(normal, (abs(normal.x) < 0.5) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0));
  vec3 tdir = cross(normal, sdir);
  return r * cos(angle) * sdir + r * sin(angle) * tdir + sqrt(1.0 - u) * normal;
}

float random(vec3 scale, float seed) {
  return fract(sin(dot(gl_FragCoord.xyz + seed, scale)) * 43758.5453 + seed);
}

vec3 radiance(Ray ray, float seed) {
  vec3 color = vec3(0.0);
  vec3 reflectance = vec3(1.0);
  for (int depth = 0; depth < 8; depth++) {
    Hit hit = intersect(ray);
    if (!(hit.distance < kInf))
      break;
    seed += float(depth);
    color += reflectance * hit.sphere.colorEmission.w;
    reflectance *= hit.sphere.colorEmission.xyz;
    vec3 hitPosition = ray.origin + ray.direction * hit.distance;
    vec3 hitNormal = (hitPosition - hit.sphere.positionRadius.xyz) / hit.sphere.positionRadius.w;
    bool bInto = dot(hitNormal, ray.direction) < 0.0;
    vec3 normal = bInto ? hitNormal : -hitNormal;
    const float kNC = 1.0;
    const float kNT = 1.5;
    vec3 reflection = reflect(ray.direction, normal);
    vec3 refraction = refract(ray.direction, normal, bInto ? (kNC / kNT) : (kNT / kNC));
    if (hit.sphere.reflectionType == kDiffuse)
      ray = Ray(hitPosition, cosineWeightedDirection(random(hitPosition, seed), random(normal, seed), normal));
    else if (any(bvec2(hit.sphere.reflectionType == kSpecular, !any(bvec3(refraction)))))
      ray = Ray(hitPosition, reflection);
    else {
      const float kR0 = (kNT - kNC) * (kNT - kNC) / ((kNT + kNC) * (kNT + kNC));
      float c = 1.0 - (bInto ? -dot(ray.direction, normal) : dot(refraction, -normal));
      float Re = kR0 + (1.0 - kR0) * c * c * c * c * c;
      float P = 0.25 + 0.5 * Re;
      if (random(refraction, seed) < P) {
        reflectance = reflectance * Re / P;
        ray = Ray(hitPosition, reflection);
      } else {
        reflectance = reflectance * (1.0 - Re) / (1.0 - P);
        ray = Ray(hitPosition + refraction * 0.01, refraction);
      }
    }
  }
  return color;
}

void main() {
  vec3 pos, ray;

  if (!setup_ray(eye, dir, pos, ray)) {
    gl_FragColor = vec4(0.);
    gl_FragDepth = 0.;
    return;
  }

  const float kGamma = 2.2;

  // Fetch current accumulation.
  vec4 current = texture2D(bg_texture, gl_FragCoord.xy / iResolution);
  vec3 curCol = pow(current.rgb, vec3(kGamma));

  float currentWeight = float(bg_weight);

  vec3 col = mix(
      radiance(Ray(pos, ray), time),
      curCol,
      currentWeight / (currentWeight + 1.));

  gl_FragColor = vec4(pow(col, vec3(1./kGamma)), 0.);
  gl_FragDepth = 0.;
}
