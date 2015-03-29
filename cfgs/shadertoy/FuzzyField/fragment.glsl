// Fuzzy Field by eiffie (https://www.shadertoy.com/view/Md2Gzh)
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// just trying to find a use for the new DoF renderer

#include "setup.inc"
#line 7

uniform vec3 par[20];
uniform int iters;
uniform int max_steps;

vec2 size = vec2(xres, yres);

//const float focalDistance=1.5,aperature=0.01,fudgeFactor=0.9;
#define aperature par[1].x //{min=.001 max=0.5 step=.001}
#define focalDistance par[1].y  //{min=0 max=5 step=.01}
#define fudgeFactor par[2].y //{min=.001 max=1 step=.001}

#define TAO 6.283
// i got these from knighty and/or darkbeam
void Rotate(inout vec2 v,float angle) {v=cos(angle)*v+sin(angle)*vec2(v.y,-v.x);}
void Kaleido(inout vec2 v,float power){Rotate(v,floor(.5+atan(v.x,-v.y)*power/TAO)*TAO/power);}
float linstep(float a, float b, float t){return clamp((t-a)/(b-a),0.,1.);}

float rand(vec2 co){// implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
  return fract(sin(dot(co*0.123,vec2(12.9898,78.233))) * 43758.5453);
}

float DE(vec3 z0)
{//a field of fuzzies
  vec4 z = vec4(z0,1.0);z.y=-z.y;//oops upsidedown
  float r=length(z.xz);
  z.xyz+=sin(vec3(z.z,length(z.xz)*0.5+(z.z+z.x)*0.5,z.x))*0.25;
  z.xz=abs(mod(z.xz,3.08)-1.54)-0.77;
  float d=max(abs(z.y)-1.0,length(z.xz+vec2(sin(z.y*3.14159)*0.04)))-0.015*clamp(1.5-abs(z.y),0.0,1.0);
  z.xyz+=sin(z.zxy*25.0)*0.01;
  z.y+=1.0;
  z*=3.0;
  Kaleido(z.xy,8.0);
  Kaleido(z.zy,8.0);
  d=min(d,(max(abs(z.y)-1.0,length(z.xz+vec2(sin(z.y*3.14159)*0.04)))-0.01*clamp(1.0-abs(z.y),0.0,1.0))*0.333);
  z.y+=1.0;
  z*=3.0;
  Kaleido(z.xy,8.0);
  Kaleido(z.zy,8.0);
  d=min(d,(max(abs(z.y)-1.0,length(z.xz+vec2(sin(z.y*3.14159)*0.04)))-0.01*clamp(1.0-abs(z.y),0.0,1.0))*0.111);
  return d;
}
vec3 mcol;
float CE(vec3 z0){mcol+=vec3(0.7)+sin(z0*10.0)*0.05;return DE(z0);}

float pixelSize;
float CircleOfConfusion(float t){//calculates the radius of the circle of confusion at length t
  return abs(focalDistance-t)*aperature+pixelSize*(1.0+t);
}
mat3 lookat(vec3 fw,vec3 up){
  fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}

float de_for_host(vec3 p) { return DE(p); }

void main() {
  pixelSize=1.0/size.y;
  vec3 ro=vec3(time*0.7+cos(time),0.5+sin(time*0.7)*0.5,time+sin(time));
  vec3 rd=lookat(vec3(1.0,1.7-ro.y*0.75,0.7),vec3(0.0,1.0,0.0))*normalize(vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,1.0));
  vec3 L=normalize(ro+vec3(0.5,2.5,0.5));
  if (!setup_ray(eye, dir, ro, rd)) {  // boxplorify view
    return;
  }
  vec4 col=vec4(0.0);//color accumulator
  float t=0.0;//distance traveled
  for(int i=1;i<48;i++){//march loop
    if(col.w>0.9 || t>20.0)continue;//bail if we hit a surface or go out of bounds
    float rCoC=CircleOfConfusion(t);//calc the radius of CoC
    float d=DE(ro)+0.5*rCoC;
    if(d<rCoC){//if we are inside add its contribution
      vec3 p=ro-rd*abs(d-rCoC);//back up to border of CoC
      mcol=vec3(0.0);//clear the color trap, collecting color samples with normal deltas
      vec2 v=vec2(rCoC*0.5,0.0);//use normal deltas based on CoC radius
      vec3 N=normalize(vec3(-CE(p-v.xyy)+CE(p+v.xyy),-CE(p-v.yxy)+CE(p+v.yxy),-CE(p-v.yyx)+CE(p+v.yyx)));
      vec3 scol=mcol*0.1666*(0.7+0.3*dot(N,L));
      scol+=pow(max(0.0,dot(reflect(rd,N),L)),8.0)*vec3(1.0,0.5,0.0);
      float alpha=fudgeFactor*(1.0-col.w)*linstep(-rCoC,rCoC,-d)/max(1.0,t);//calculate the mix like cloud density
      col+=vec4(scol*alpha,alpha);//blend in the new color
    }
    d=abs(fudgeFactor*d*(0.8+0.2*rand(gl_FragCoord.xy*vec2(i))));//add in noise to reduce banding and create fuzz
    ro+=d*rd;//march
    t+=d;
  }//mix in background color
  vec3 scol=mix(vec3(0.025,0.1,0.05)+rd*0.025,vec3(0.1,0.2,0.3)+rd*0.1,smoothstep(-0.1,0.1,rd.y));
  col.rgb+=scol*(1.0-clamp(col.w,0.0,1.0));

  //gl_FragColor = vec4(clamp(col.rgb,0.0,1.0),1.0);
  write_pixel(dir, t, clamp(col.rgb, 0.0, 1.0));  // boxplorify write
}
