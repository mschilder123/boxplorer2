//by @paulofalcao, https://www.shadertoy.com/view/Xdj3Dt

#include "setup.inc"
#line 5
float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);

float t=iGlobalTime;
vec2 iRes=iResolution.xy;

float stime=sin(t);
float ctime=cos(t);

float inObj(in vec3 p){
  float oP=length(p);
  p.x=sin(p.x)+stime;
  p.z=sin(p.z)+ctime;
  return float(min(length(p)-1.5-sin(oP-t*4.0),p.y+3.0));
}

void main(void){
  vec2 iRes=vec2(iRes);
  vec2 vPos=-1.0+2.0*gl_FragCoord.xy/iRes;

  float zf=cos(t*0.2)*0.4+0.6;
  vec3 vuv=vec3(stime,1,0);
  vec3 vrp=vec3(sin(t*0.7)*10.0,0,cos(t*0.9)*10.0)*zf; 
  vec3 prp=vec3((sin(t*0.7)*20.0+20.0)*zf,stime*4.0+4.0+3.0,(cos(t*0.6)*20.0+14.0)*zf)+vrp;

  vec3 vpn=normalize(vrp-prp);
  vec3 u=normalize(cross(vuv,vpn));
  vec3 v=cross(vpn,u);
  vec3 scrCoord=vpn+vPos.x*u*iRes.x/iRes.y+vPos.y*v;
  vec3 scp=normalize(scrCoord);

  if (!setup_ray(eye, dir, prp, scp)) return;

  const vec3 e = vec3(0.1,0,0);
  const float maxd=80.0;

  float s=0.1;
  vec3 c,p,n;

#define FL -2.5
#define CL 3.5
  float f=-(prp.y-CL)/scp.y;
  if (f>0.0)
    p=prp+scp*f;
  else f=maxd;

  vec3 outc=vec3(0,0,0);

  float far=0.0;
  for (int ref=0;ref<=1;ref++){
    if (ref>0){
      scp=reflect(scp,n);
      prp=p;
      s=0.1;f=0.1;
    }
    for(int i=0;i<32;i++){
      f+=s;
      p=prp+scp*f;
      s=inObj(p);
      if (abs(s)<.01||f>maxd||p.y>CL||(ref>0&&i>16)) break;
    }
   
    if (f<maxd&&p.y<CL){
      if(p.y<FL){  // draw tile
        if (fract(p.x/4.0)>.5)
          if (fract(p.z/4.0)>.5)
            c=vec3(0,0,0);
          else
            c=vec3(1,1,1);
        else
          if (fract(p.z/4.0)>.5)
            c=vec3(1,1,1);
          else
            c=vec3(0,0,0);
        c=c*max(inObj(vec3(p.x,p.y+1.0,p.z)),0.5);
        n=vec3(0,1,0);
      } else {  // draw blob
        float d=length(p);
        c=vec3((sin(d*.25-t*4.0)+1.0)/2.0,
               (stime+1.0)/2.0,
               (sin(d-t*4.0)+1.0)/2.0);
        n=normalize(
          vec3(s-inObj(p-e.xyy),
               s-inObj(p-e.yxy),
               s-inObj(p-e.yyx)));
      }
      float b=dot(n,normalize(prp-p));
      if (ref==0) {
        outc=((b+0.2)*c+pow(b,54.0))*0.7;
        far=1.0-f*.01;
	    } else {
		    if (prp.y>FL+.1) outc+=(b+0.2)*c*0.3;
	    }
    }
    else break;
  }
  write_pixel(dir, far, outc*far);
}
