//Simple path tracer by Movania Muhammad Mobeen
//Hacked up by marius to work /w boxplorer2

uniform float xres, yres, speed, time;
varying vec3 eye, dir;  
uniform sampler2D bg_texture;
uniform vec3 par[20];

#include "setup.inc"

vec2 TEX_SIZE = vec2(xres, yres);

vec3 light = par[0];

struct Sphere {
  vec4 center_radius;
  vec3 surfaceColor;
};

struct Box {
  vec3 min, max;
};

Box rootCube = Box(par[1], par[2]);

Sphere sphere0 = Sphere(vec4(par[3], .12), vec3(1.,.1,.1));
Sphere sphere1 = Sphere(vec4(par[4], .12), vec3(.9));
Sphere sphere2 = Sphere(vec4(par[5], .12), vec3(.9));
Sphere sphere3 = Sphere(vec4(par[6], .12), vec3(.1,1.,.1));
 
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
  if(hit.x < cube.min.x + 0.0001) return vec3(-1.0, 0.0, 0.0);  
  else if(hit.x > cube.max.x - 0.0001) return vec3(1.0, 0.0, 0.0);  
  else if(hit.y < cube.min.y + 0.0001) return vec3(0.0, -1.0, 0.0); 
  else if(hit.y > cube.max.y - 0.0001) return vec3(0.0, 1.0, 0.0);    
  else if(hit.z < cube.min.z + 0.0001) return vec3(0.0, 0.0, -1.0); 
  else return vec3(0.0, 0.0, 1.0);  
}

float intersectSphere(vec3 origin, vec3 ray, Sphere s) {  
  vec3 toSphere = origin - s.center_radius.xyz;
  float sphereRadius = s.center_radius.w;
  float a = dot(ray, ray);
  float b = 2.0 * dot(toSphere, ray);
  float c = dot(toSphere, toSphere) - sphereRadius*sphereRadius;
  float discriminant = b*b - 4.0*a*c;
  if(discriminant >= 0.0) {
    float t = (-b - sqrt(discriminant)) / (2.0 * a);
    if(t >= 0.0) return t;
  }
  return 10000.0;
}

vec3 normalForSphere(vec3 hit, Sphere s) {
  return normalize(hit - s.center_radius.xyz);
}

// [0..1>
float random(vec3 scale, float seed) {
  return fract(sin(dot(gl_FragCoord.xyz + seed, scale)) * 43758.5453 + seed);
}

// normalized
vec3 uniformlyRandomDirection(float seed) {
  float u = random(vec3(12.9898, 78.233, 151.7182), seed);
  float v = random(vec3(63.7264, 10.873, 623.6736), seed);
  float z = 1.0 - 2.0 * u;
  float r = sqrt(1.0 - z * z);
  float angle = 6.283185307179586 * v;
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
    vec2 tRoom = intersectCube(origin, ray, rootCube);       
    float t = 10000.0;          

    if(tRoom.x < tRoom.y) {
      t = tRoom.y;  
    } else {
      // missed the room
      break;
    }

    float tSphere0 = intersectSphere(origin, ray, sphere0);
    float tSphere1 = intersectSphere(origin, ray, sphere1);
    float tSphere2 = intersectSphere(origin, ray, sphere2);
    float tSphere3 = intersectSphere(origin, ray, sphere3);

    Sphere sphere;
    if(tSphere0 < t) { t = tSphere0;sphere=sphere0;}
    if(tSphere1 < t) { t = tSphere1;sphere=sphere1;} 
    if(tSphere2 < t) { t = tSphere2;sphere=sphere2;} 
    if(tSphere3 < t) { t = tSphere3;sphere=sphere3;} 

    if (t == 1000.0) break;

    vec3 hit = origin + ray * t;
    vec3 surfaceColor = vec3(0.75);
    vec3 normal;

    if(t == tRoom.y) {
      normal = -normalForCube(hit, rootCube); 
      if(hit.x < -0.9999) {
        surfaceColor = vec3(0.1, 0.4, 1.0); 
      } else if(hit.x > 0.9999) {
        surfaceColor = vec3(1.0, 0.9, 0.1); 
      }
    } else {
      normal = normalForSphere(hit, sphere);
      surfaceColor = sphere.surfaceColor;
    }

    // pick scatter ray, no bias
    vec3 newray = uniformlyRandomDirection(time + float(bounce));
    if(dot(normal, newray) < 0.0) newray = -newray;      

    vec3 toLight = light - hit;     
    vec3 ldir = normalize(toLight);

    // blinn-phong
    float diffuse = max(0.0, dot(ldir, normal));
    vec3 halfLV = normalize(ldir + -ray);
    float specularHighlight = pow(max(dot(normal, halfLV), 0.0), 32.0);

    float shadowIntensity = shadow(hit + normal * 0.001, toLight);

    colorMask *= surfaceColor;
    accumulatedColor += colorMask * (diffuse * shadowIntensity);      
    accumulatedColor += specularHighlight * shadowIntensity;

    origin = hit + newray * .001;  // step away a bit from hit
    ray = newray;
  }

  return accumulatedColor * 0.5;  
}

void main() {
   vec3 pos, ray;

   if (!setup_ray(eye, dir, pos, ray)) {
     gl_FragColor = vec4(0.);
     gl_FragDepth = 0.;
     return;
   }

  // Jitter light a bit
  vec3 newLight = light + uniformlyRandomVector(time + 53.0) * 0.1;   

  // Fetch current accumulation
  vec4 current = texture2D(bg_texture, gl_FragCoord.xy / TEX_SIZE);

  float textureWeight = 999.*current.w;  // [0..999]

  vec3 col = mix(
      calculateColor(pos, ray, newLight),
      current.rgb,
      textureWeight / (textureWeight + 1.));

  gl_FragColor = vec4(col, (textureWeight + 1.) / 1000.);
  gl_FragDepth = 0.;
}
