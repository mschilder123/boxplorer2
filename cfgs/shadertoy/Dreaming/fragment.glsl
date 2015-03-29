// By Kali, from https://www.shadertoy.com/view/ldXGDl

#include "setup.inc"
#line 5

float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);

void main(void)
{
	vec2 uv = gl_FragCoord.xy / iResolution.xy - .5;
	uv.x*=iResolution.x/iResolution.y;
	//if (iMouse.z>0.) uv+=(iMouse.xy/iResolution.xy-.5)*2.;
	float t=iGlobalTime*.2;
	float fov=1.5+sin(t);
	vec3 ddir=normalize(vec3(uv*fov*mat2(cos(t),sin(t),-sin(t),cos(t)),1.));
	vec3 from=vec3(t*-5.,sin(t*1.28352),0.)*.005;
  if (!setup_ray(eye, dir, from, ddir)) {  // boxplorify view
    return;
  }
	float dist=-0.015; 
	vec3 vol=vec3(0.);
	for (int v=0; v<500; v++) {
		dist+=.00015;
		vec3 p=from+dist*ddir*vec3(vec2(sign(dist)),1.);				
		vec3 disp=texture2D(iChannel0,vec2(length(p.xy),p.z*.3)*3.2568).xyz-.5;
		vol+=pow(length(abs(.5-mod(p+disp*.013,vec3(.01))/.01)),14.)
			*(.7+normalize(disp)*.3)*exp(-500.*dist*dist);
	}
	vol=(vol+.2)*min(1.,t-.05)*vec3(1.,.8,.7);
  write_pixel(dir, 1.0, pow(vol, vec3(1.3)));  // boxplorify write
}
