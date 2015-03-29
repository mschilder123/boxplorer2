//After http://madebyevan.com/webgl-path-tracing/webgl-path-tracing.js
//Hacked up by marius to work /w boxplorer2

#include "setup.inc"
#line 6

uniform vec3 par[20];

vec2 iResolution = vec2(xres, yres);

vec3 light = par[0];

struct Sphere {
  vec4 center_radius;
  vec3 surfaceColor;
  float glossiness;
  float specular;
};

struct Box {
  vec3 min, max;
};

const float infinity = 10000.0;
const float epsilon = 0.0001;
const float TWOPI = 6.283185307179586;
const float wallGlossiness = .7;

Box rootCube = Box(par[1], par[2]);

Sphere sphere0 = Sphere(vec4(par[3], .25), vec3(1.,.1,.1), .5, 3.);
Sphere sphere1 = Sphere(vec4(par[4], .25), vec3(.9), .0, 20.);
Sphere sphere2 = Sphere(vec4(par[5], .25), vec3(.9), .0, 20.);
Sphere sphere3 = Sphere(vec4(par[6], .25), vec3(.1,1.,.1), .5, 3.);

vec2 intersectCube(vec3 origin, vec3 ray, Box cube) {   
  vec3   tMin = (cube.min - origin) / ray;    
  vec3   tMax = (cube.max - origin) / ray;    
  vec3     t1 = min(tMin, tMax);    
  vec3     t2 = max(tMin, tMax);
  float tNear = max(max(t1.x, t1.y), t1.z);
  float  tFar = min(min(t2.x, t2.y), t2.z);
  return vec2(tNear, tFar); 
}

vec3 normalForCube(vec3 hit, Box cube)  { 
  if(hit.x < cube.min.x + epsilon) return vec3(-1.0, 0.0, 0.0);  
  else if(hit.x > cube.max.x - epsilon) return vec3(1.0, 0.0, 0.0);  
  else if(hit.y < cube.min.y + epsilon) return vec3(0.0, -1.0, 0.0); 
  else if(hit.y > cube.max.y - epsilon) return vec3(0.0, 1.0, 0.0);    
  else if(hit.z < cube.min.z + epsilon) return vec3(0.0, 0.0, -1.0); 
  else return vec3(0.0, 0.0, 1.0);  
}

bool insideSphere(vec3 origin, Sphere s) {
  float toSphere = length(origin - s.center_radius.xyz);
  return toSphere < s.center_radius.w;
}

float intersectSphere(vec3 origin, vec3 ray, Sphere s) {  
  vec3 toSphere = origin - s.center_radius.xyz;
  float sphereRadius = s.center_radius.w;
  float a = dot(ray, ray);
  float b = 2.0 * dot(toSphere, ray);
  float c = dot(toSphere, toSphere) - sphereRadius*sphereRadius;
  float discriminant = b*b - 4.0*a*c;
  if(discriminant > 0.0) {
    float t = (-b - sqrt(discriminant)) / (2.0 * a);
    if(t > 0.0) return t;
  }
  return infinity;
}

vec3 normalForSphere(vec3 hit, Sphere s) {
  return normalize(hit - s.center_radius.xyz);
}

// [0..1>
float random(vec3 scale, float seed) {
  return fract(sin(dot(gl_FragCoord.xyz + seed, scale)) * 43758.5453 + seed);
}

vec3 cosineWeightedDirection(float seed, vec3 normal) {
  float u = random(vec3(12.9898, 78.233, 151.7182), seed);
  float v = random(vec3(63.7264, 10.873, 623.6736), seed);
  float r = sqrt(u);
  float angle = TWOPI * v;
  // compute basis from normal
  vec3 sdir, tdir;
  if (abs(normal.x)<.5) {
    sdir = cross(normal, vec3(1,0,0));
  } else {
    sdir = cross(normal, vec3(0,1,0));
  }
  tdir = cross(normal, sdir);
  return r*cos(angle)*sdir + r*sin(angle)*tdir + sqrt(1.-u)*normal;
}

// normalized
vec3 uniformlyRandomDirection(float seed) {
  float u = random(vec3(12.9898, 78.233, 151.7182), seed);
  float v = random(vec3(63.7264, 10.873, 623.6736), seed);
  float z = 1.0 - 2.0 * u;
  float r = sqrt(1.0 - z * z);
  float angle = TWOPI * v;
  return vec3(r * cos(angle), r * sin(angle), z);
}

// length |[0..1>|
vec3 uniformlyRandomVector(float seed) {
  return uniformlyRandomDirection(seed) *
         (random(vec3(36.7539, 50.3658, 306.2759), seed));
}

float shadow(vec3 origin, vec3 ray) {
  float tSphere0 = intersectSphere(origin, ray, sphere0);
  if(tSphere0 < 1.0) return 0.0;
  float tSphere1 = intersectSphere(origin, ray, sphere1);
  if(tSphere1 < 1.0) return 0.0;
  float tSphere2 = intersectSphere(origin, ray, sphere2);
  if(tSphere2 < 1.0) return 0.0;
  float tSphere3 = intersectSphere(origin, ray, sphere3);
  if(tSphere3 < 1.0) return 0.0;
  return 1.0;
}

vec3 calculateColor(vec3 origin, vec3 ray, vec3 light) {    
  vec3 colorMask = vec3(1.0);   
  vec3 accumulatedColor = vec3(0.0);

  for(int bounce = 0; bounce < 5; bounce++) {     
    float t = infinity;          

    vec2 tRoom = intersectCube(origin, ray, rootCube);       
    if(tRoom.x < tRoom.y) t = tRoom.y;  

    float tSphere0 = intersectSphere(origin, ray, sphere0);
    float tSphere1 = intersectSphere(origin, ray, sphere1);
    float tSphere2 = intersectSphere(origin, ray, sphere2);
    float tSphere3 = intersectSphere(origin, ray, sphere3);

    Sphere sphere;
    if(tSphere0 < t) { t = tSphere0;sphere=sphere0;}
    if(tSphere1 < t) { t = tSphere1;sphere=sphere1;} 
    if(tSphere2 < t) { t = tSphere2;sphere=sphere2;} 
    if(tSphere3 < t) { t = tSphere3;sphere=sphere3;} 

    vec3 hit = origin + ray * t;
    vec3 toLight = light - hit;
    vec3 surfaceColor = vec3(0.75);
    float specularHighlight = 0.;
    vec3 normal;
    vec3 newray;

    if(t == tRoom.y) {
      normal = -normalForCube(hit, rootCube); 
      if(hit.x < -0.9999) {
        surfaceColor = vec3(0.1, 0.4, 1.0);  // blue
        // diffuse ray
        newray = cosineWeightedDirection(time + float(bounce), normal);
      } else if(hit.x > 0.9999) {
        surfaceColor = vec3(1.0, 0.9, 0.1);  // yellow
        // diffuse ray
        newray = cosineWeightedDirection(time + float(bounce), normal);
      } else {
        // glossy ray
        newray = normalize(reflect(ray, normal)) +
                      uniformlyRandomVector(time + float(bounce)) * wallGlossiness;
        vec3 reflectedLight = normalize(reflect(toLight, normal));
        specularHighlight = max(0.0,
                                dot(reflectedLight, normalize(hit - origin)));
        specularHighlight = pow(specularHighlight, 3.0);
      }
    } else {
      normal = normalForSphere(hit, sphere);
      surfaceColor = sphere.surfaceColor;
      // glossy / reflective ray
      newray = normalize(
          reflect(ray, normal)) + uniformlyRandomVector(time + float(bounce))
                      * sphere.glossiness;
      vec3 reflectedLight = normalize(reflect(toLight, normal));
      specularHighlight = max(0.0,
                              dot(reflectedLight, normalize(hit - origin)));
      specularHighlight = pow(specularHighlight, sphere.specular);
    }

    float diffuse = max(0.0, dot(normalize(toLight), normal));
    float shadowIntensity = shadow(hit + normal * epsilon, toLight);

    colorMask *= surfaceColor;
    accumulatedColor += colorMask * (.5 * diffuse * shadowIntensity);
    accumulatedColor += colorMask * specularHighlight * shadowIntensity;

    origin = hit + newray * epsilon;  // step away a bit from hit
    ray = newray;
  }

  return accumulatedColor;
}

void main() {
  vec3 pos, ray;

  if (!setup_ray(eye, dir, pos, ray)) {
    gl_FragColor = vec4(0.);
    gl_FragDepth = 0.;
    return;
  }

  const float kGamma = 1.0;

  // Fetch current accumulation
  vec4 current = texture2D(iBackbuffer, gl_FragCoord.xy / iResolution);
  vec3 curCol = pow(current.xyz, vec3(kGamma));

  // Jitter light a bit for soft shadows.
  vec3 newLight = light + uniformlyRandomVector(time + 53.0) * 0.3;

  float currentWeight = float(iBackbufferCount);

  vec3 col = mix(
      calculateColor(pos, ray, newLight),
      curCol,
      currentWeight / (currentWeight + 1.));

  gl_FragColor = vec4(pow(col, vec3(1./kGamma)), 0.);
  gl_FragDepth = 0.;
}
