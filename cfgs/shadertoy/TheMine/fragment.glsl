// created by Vinicius Graciano Santos - vgs/2014
// http://vgsx.wordpress.com/2014/07/08/making-of-the-mine-modeling/

#include "setup.inc"
#line 6

#define iChannel1 iChannel0
#define iChannel2 iChannel0
#define iChannel3 iChannel0

float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);

// increase/decrease the STEPS value if you have a fast/slow gpu
#define STEPS 100
#define EPS 0.01

// uncomment to enable bump mapping (crashes on ANGLE-based browsers)
//#define BUMP 5.0

// iq's 3D noise function: https://www.shadertoy.com/view/4sfGzS.
float hash(float n ) { return fract(sin(n)*43758.5453123); }
float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    float n = p.x + p.y*157.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+157.0), hash(n+158.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+270.0), hash(n+271.0),f.x),f.y),f.z);
}

float fbm(vec3 p) {
    float k = 0.0;
        
    k += 1.000*noise(p); p*=2.0;
    k += 0.500*noise(p); p*=2.0;
    k += 0.250*noise(p);
    return k/1.75;
}

vec2 track(float z) {
    // play with these constants for some fun!
    float x = cos(0.3*z);
    float y = -cos(0.2*z) - 0.2*sin(0.8*z - 2.0);
    return vec2(x, y);
}

float cave(vec3 p) {
    const float k = 4.0;
    return 1.6-pow(pow(abs(p.x), k) + pow(abs(p.y), k), 1.0/k);
}

float box(vec2 p, vec2 b, float r) {
    return length(max(abs(p)-b, 0.0)) - r;
}

float box(vec3 p, vec3 b, float r) {
    return length(max(abs(p)-b, 0.0)) - r;
}

float support(vec3 p) {
    const vec4 c = vec4(0.15, 0.2, 2.0, 1.2);
    vec3 q = vec3(abs(p.x) - c.z, p.y - c.z, mod(p.z, 6.0) - 3.0);
    float d = box(q.xz, c.xx, 0.05);
    d = min(d, box(q.yz, c.yx, 0.05)); 
    q.x += q.y + c.w;
    return min(d, box(q.xz, c.xx, 0.05));
}

float plank(vec3 p) {
    vec3 q = vec3(p.x, p.y + 1.9, mod(p.z, 2.0) - 1.0);
    return box(q, vec3(1.5, 0.05, 0.2), 0.01);
}

float rails(vec3 p) {
    vec2 q = vec2(abs(p.x)-1.0, p.y + 1.7);
    float d = box(q, vec2(0.1), 0.01); q.x += 0.2;
    return max(d, 0.11 - length(q));
}

vec2 dist_field(vec3 p) {
    p.xy += track(p.z);
    vec2 res = vec2(cave(p) + fbm(p), 0.0);
    
    float d = support(p);
    if (d < res.x) res = vec2(d, 1.0);
    d = plank(p);
    if (d < res.x) res = vec2(d, 1.0);
    d = rails(p);
    if (d < res.x) res = vec2(d, 2.0);
    return res;
}

vec3 normal(vec3 p) {
    vec2 q = vec2(0.01, 0.0);
    return normalize(vec3(dist_field(p+q.xyy).x - dist_field(p-q.xyy).x,
                          dist_field(p+q.yxy).x - dist_field(p-q.yxy).x,
                          dist_field(p+q.yyx).x - dist_field(p-q.yyx).x));
}

vec3 cubeMap(sampler2D samp, vec3 q, vec3 n) {
    vec3 x = texture2D(samp, q.zy).rgb;
    vec3 y = texture2D(samp, q.zx).rgb;
    vec3 z = texture2D(samp, q.xy).rgb;
    return abs(n.x)*x + abs(n.y)*y + abs(n.z)*z;
}

#ifdef BUMP
vec3 normalMap(sampler2D samp, vec2 q) {
    vec2 p = vec2(0.01, 0.0);
    vec3 a = BUMP*(texture2D(samp, q+p.xy).rgb - texture2D(samp, q-p.xy).rgb);
    vec3 b = BUMP*(texture2D(samp, q+p.yx).rgb - texture2D(samp, q-p.yx).rgb);
    return normalize(cross(vec3(1.0, 0.0,  (a.x+a.y+a.z)/3.0),
                           vec3(0.0, 1.0, (b.x+b.y+b.z)/3.0)));
}

vec3 cubeNormalMap(sampler2D samp, vec3 q, vec3 n) {
    vec3 x = normalMap(samp, q.zy);
    vec3 y = normalMap(samp, q.zx);
    vec3 z = normalMap(samp, q.xy);
    return abs(n.x)*x + abs(n.y)*y + abs(n.z)*z;
}
#endif

vec3 shade(in vec3 ro, in vec3 rd, float t, float id) {
    vec3 key_l = -rd;
    vec3 key_c = vec3(243.0, 141.0, 21.0)/25.5;
    
    vec3 fill_l = vec3(0.0, 0.0, -1.0);
    vec3 fill_c = 0.2*key_c;
    
    vec3 q = ro + t*rd;
    vec3 n = normal(q);

    vec3 mat = vec3(1.0); float shin = 0.0;
    if (id == 0.0) {
        shin = 25.0;
        mat = mix(cubeMap(iChannel0, q, n), cubeMap(iChannel1, q, n), noise(q));
    } else if (id == 1.0) {
        shin = 50.0;
        mat = cubeMap(iChannel2, q, n).ggg;
        #ifdef BUMP
        n = normalize(n + cubeNormalMap(iChannel2, q, n));
        #endif
    } else {
        shin = 75.0;
        mat = cubeMap(iChannel3, q, n).bbb;
        #ifdef BUMP
        n = normalize(n + cubeNormalMap(iChannel3, q, n));
        #endif
    }
    mat = pow(abs(mat), vec3(2.2));
    vec3 col = vec3(0.0);
    
    // key light.
    float lamb = max(0.0, dot(n, key_l));
    float spec = lamb > 0.0 ? pow(max(0.0, dot(n, normalize(key_l-rd))), shin) : 0.0;
    col += key_c*mat*(0.6*lamb + 0.4*spec)*pow(max(-dot(-rd,normalize(q)), 0.0), 10.0)/(0.1*t*t);
    
    // fill light.
    lamb = max(0.0, dot(n, fill_l));
    spec = lamb > 0.0 ? pow(max(0.0, dot(n, normalize(fill_l-rd))), shin) : 0.0;
    col += fill_c*mat*(0.6*lamb + 0.4*spec)/(0.4*t*t);
    
    col = mix(vec3(0.05), col, exp(-0.05*t));
    return col/(col+1.0);
    
}

void main(void) {
    vec2 uv = -1.0 + 2.0*gl_FragCoord.xy / iResolution.xy;
    uv.x *= iResolution.x/iResolution.y;
    
    vec3 ro = vec3(0.0, 0.0, -10.0*iGlobalTime-10.0);
    ro.xy -= track(ro.z);
    vec3 rd = normalize(vec3(uv, -1.0));

    vec3 dummy_ro;
    if (!setup_ray(eye, dir, dummy_ro, rd)) return;
    ro += (eye - dummy_ro);
    
    float t = 0.0; vec2 res = vec2(3.0);
    for (int i = 0; i < STEPS; ++i) {
        if (res.x < EPS || t > 32.0) continue;
        res = dist_field(ro + t*rd); t += 0.7*res.x;
    }
    
    vec3 col = shade(ro, rd, t, res.y);
    col = pow(abs(col), vec3(1.0/2.2));
    col = smoothstep(0.0, 1.0, col);

    write_pixel(dir, t, col);
}
