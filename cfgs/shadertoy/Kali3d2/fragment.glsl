// From https://www.shadertoy.com/view/Xs2GDK

// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

//Based on https://www.shadertoy.com/view/4ds3zn by IQ

#include "setup.inc"
#line 8

uniform int iters, max_steps;
uniform vec3 par[20];

vec2 iResolution = vec2(xres, yres);
vec3 iMouse = vec3(0., 0., 0.);
float iGlobalTime = time;

vec4 ot; 

//vec3 C =  //vec3(1.05);
 //vec3(.7,.9,1.41);
#define CVector  par[2]
float igt = iGlobalTime;
float kzoom=1.2;

float map( vec3 p )
{
  float dr = 1.0;
  
  ot = vec4(1000.0); 
  float r2;
  vec3 C = (abs(CVector));
  
  for( int i=0; i<iters; i++ )
  {
    r2 = dot(p, p);
    if(r2>100.) break;

    ot = min( ot, vec4(abs(p),r2) );

    //Kali formula 
    p= abs(p)/r2-C; 

    dr= dr/r2;  // _not_ /length(p)
  } 
  return .1*(abs(p.x)+abs(p.y))*length(p)/dr;
  //return .15*length(p.xz)*length(p.xy)/dr;
  //return .125*sqrt(r2)*log(r2)/dr;
  //return .1*length(p)/dr;
}

//float de_for_host(vec3 p) { return map(p); }
const float maxd = 20.;

float trace( in vec3 ro, in vec3 rd )
{
  float precis = 0.001;
      
  float h=precis*2.0;
  float t = 0.0;
  for( int i=0; i<max_steps; i++ )
  {
    if( t>maxd || h<precis*(.1+t)) break;

    t += h;
    h = map( ro+rd*t );
  }

  //if( t>maxd ) t=-1.0;
  return t;
}

vec3 calcNormal( in vec3 pos )
{
  vec3  eps = vec3(.0001,0.0,0.0);
  vec3 nor;
  nor.x = map(pos+eps.xyy) - map(pos-eps.xyy);
  nor.y = map(pos+eps.yxy) - map(pos-eps.yxy);
  nor.z = map(pos+eps.yyx) - map(pos-eps.yyx);
  return normalize(nor);
}

void main(void)
{
  vec2 p = -1.0 + 2.0*gl_FragCoord.xy / iResolution.xy;
  p.x *= iResolution.x/iResolution.y;

  vec2 m = vec2(-0.5)*6.28;
  if( iMouse.z>0.0 )m = (iMouse.xy/iResolution.xy-.5)*6.28;
  m+=.5*vec2(cos(0.15*igt),cos(0.09*igt))+.3;      
  
    // camera

  vec3 ta = vec3(0.,.2*sin(0.12*igt),0.);
  vec3 ro = ta- kzoom*vec3( cos(m.x)*cos(m.y), sin(m.y), sin(m.x)*cos(m.y));
  
  vec3 cw = normalize(ta-ro);
  vec3 cp = vec3(0.,1.,0.0);
  vec3 cu = normalize(cross(cw,cp));
  vec3 cv = normalize(cross(cu,cw));
  vec3 rd = normalize( p.x*cu + p.y*cv + 2.0*cw );

  if (!setup_ray(eye, dir, ro, rd)) {  // boxplorify view
    return;
  }

  // trace  
  vec3 col = vec3(0.8,0.8,1.);
  float t = trace( ro, rd );
  if( t<=maxd ) {
    vec3 pos = ro + t*rd;
    vec3 nor = calcNormal( pos );
    
    // lighting
    vec3 light1 = vec3(  0.577, 0.577, -0.577 );
    vec3 light2 = vec3( -0.707, -0.707,0.0  );
    float key = clamp( dot( light1, nor ), 0.0, 1.0 );
    float bac = clamp( 0.2 + 0.8*dot( light2, nor ), 0.0, 1.0 );
    float amb = (0.7+0.3*nor.y);
    float ao = pow( clamp(ot.w*2.0,0.2,1.0), 1.2 );   
    vec3 brdf = vec3(ao)*(.4*amb+key+.2*bac);

    // material   
    vec3 rgb = vec3(1.0);
    
    rgb =(0.4*abs(sin(2.5+(vec3(.5*ot.w,ot.y*ot.y,2.-5.*ot.w))))+1.6*sin(vec3(-0.2,-0.6,0.8)+0.+ot.x*18.))*.85 + .15;
    rgb.gbr=mix(rgb,rgb.bgr+vec3(0.3,0.1,-.2),0.5+.5*sin(8.5*ot.w));

    // color
    col = mix(vec3(0.8,0.8,1.),rgb*brdf,exp(-0.08*t));
  }

  write_pixel(dir, t, col);  // boxplorify write
}
