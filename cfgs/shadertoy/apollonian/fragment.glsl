// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//
// I can't recall where I learnt about this fractal.
//
// Coloring and fake occlusions are done by orbit trapping, as usual (unless somebody has invented
// something new in the last 4 years that i'm unaware of, that is)
//
// From: https://www.shadertoy.com/view/4ds3zn

// Camera position, direction and eps from vertex shader.
varying vec3 eye;
varying vec3 dir;
varying float zoom;

uniform float xres, yres, time, speed;
uniform int iters;    // {min=1 max=1000} Number of fractal iterations.
uniform int color_iters;    // {min=1 max=1000} Number of fractal iterations.

vec2      iResolution = vec2(xres, yres);     // viewport resolution (in pixels)
float     iGlobalTime = time;     // shader playback time (in seconds)

uniform vec3 par[1];
#define ss par[0].x  // {min=.1 max=2.0 step=.001}

#include "setup.inc"
#line 27

float de( vec3 p ) {
  float scale = 1.0;  
  for( int i=0; i<iters; i++ ) {
    p = -1.0 + 2.0*fract(0.5*p+0.5);
    float r2 = dot(p,p);      
    float k = max(ss/r2,0.1);
    p     *= k;
    scale *= k;
  }
  return 0.25*abs(p.y)/scale;
}

float de_for_host(vec3 p) { return de(p); }

vec4 trap( vec3 p ) {
  vec4 orb = vec4(1000.0);  
  for( int i=0; i<color_iters; i++ ) {
    p = -1.0 + 2.0*fract(0.5*p+0.5);
    float r2 = dot(p,p);    
    orb = min( orb, vec4(abs(p),r2) );  
    float k = max(ss/r2,0.1);
    p     *= k;
  }
  return orb;
}

float trace( in vec3 ro, in vec3 rd ) {
  float maxd = 100.0;
  float m_zoom = zoom * .5 / xres;
  float precis = 0.0001;
    float h = 0.0;
  float t = 0.0;
    for( int i=0; i<200; i++ ) {
      h = de( ro+rd*t );
    if( abs(h)<precis||t>maxd ) break;
        t += h;
      precis = m_zoom * t;
    }

    if( t>maxd ) t=-1.0;
    return t;
}

vec3 calcNormal( in vec3 pos ) {
  vec3 eps = vec3(.0001,0.0,0.0);
  vec3 nor;
  nor.x = de(pos+eps.xyy) - de(pos-eps.xyy);
  nor.y = de(pos+eps.yxy) - de(pos-eps.yxy);
  nor.z = de(pos+eps.yyx) - de(pos-eps.yyx);
  return normalize(nor);
}

void main(void) {
  vec3 ro;
  vec3 rd;
  if (!setup_ray(eye, dir, ro, rd)) {  // boxplorify view
    return;
  }

    // trace  
  vec3 col = vec3(0.0);
  float t = trace( ro, rd );
  if( t>0.0 ) {
    vec3 pos = ro + t*rd;
    vec4 tra = trap( pos );
    vec3 nor = calcNormal( pos );
    
    // lighting
    vec3  light1 = vec3(  0.577, 0.577, -0.577 );
    vec3  light2 = vec3( -0.707, 0.000,  0.707 );
    float key = clamp( dot( light1, nor ), 0.0, 1.0 );
    float bac = clamp( 0.2 + 0.8*dot( light2, nor ), 0.0, 1.0 );
    float amb = (0.7+0.3*nor.y);
    float ao = pow( clamp(tra.w*2.0,0.0,1.0), 1.2 );

    vec3 brdf  = 1.0*vec3(0.40,0.40,0.40)*amb*ao;
    brdf += 1.0*vec3(1.00,1.00,1.00)*key*ao;
    brdf += 1.0*vec3(0.40,0.40,0.40)*bac*ao;

    // material   
    vec3 rgb = vec3(tra.yzx);
    rgb = mix( rgb, vec3(1.0,0.80,0.2), clamp(6.0*tra.x,0.0,1.0) );
    rgb = mix( rgb, vec3(1.0,0.55,0.0), pow(clamp(1.0-2.0*tra.y,0.0,1.0),8.0) );

    // color
    col = rgb*brdf*exp(-0.2*t);
    }

    col = sqrt(col);
  
    col = mix( col, smoothstep( 0.0, 1.0, col ), 0.25 );
  
    write_pixel(dir, t, col);  // boxplorify write
}
