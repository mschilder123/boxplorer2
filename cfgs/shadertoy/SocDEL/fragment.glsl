// SoC with DEL by eiffie (adding Distance Estimated Light to the Sphere of Confusion renderer)
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// From https://www.shadertoy.com/view/lsjGzW

uniform float xres, yres, speed, time;
varying vec3 eye, dir;
uniform vec3 par[20];
uniform int iters;
uniform int max_steps;

vec2 size = vec2(xres, yres);

float DE(vec3);
float de_for_host(vec3 p) { return DE(p); }
#include "setup.inc"
#line 17

#define aperture par[1].x //{min=.001 max=0.5 step=.001}
#define focalDistance par[1].y  //{min=0 max=5 step=.01}
#define fudgeFactor par[2].y //{min=.001 max=1 step=.001}
//float pixelSize,focalDistance,aperture,fudgeFactor=0.78,shadowCone=0.5;
float pixelSize,shadowCone=0.5;

bool bColoring=false;
vec3 mcol;
mat2 rmx;
const vec4 p0=vec4(0.0,0.0,4.0,1.0);
const vec3 rc=vec3(2.633,0.033,2.133);
float DE(in vec3 z0){//amazing box by tglad
	vec4 z = vec4(z0,1.0);
	for (int n = 0; n < 4; n++) {
		z.xyz=clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz;
		z*=2.0/clamp(dot(z.xyz,z.xyz),0.1,1.0);
		z+=p0;
		z.xy=z.xy*rmx;
	}
	if(bColoring)mcol+=vec3(0.7,0.6,0.4)+z.xyz*0.2;
	z.xyz=max(abs(z.xyz)-rc,vec3(0.0));
	return (length(z.xyz)-0.1)/z.w;
}
float DEL(in vec3 z0){//distance to nearest light
	vec4 z = vec4(z0,1.0);
	z.xyz=clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz;
	z*=2.0/clamp(dot(z.xyz,z.xyz),0.1,1.0);
	z+=p0;
	z.xy=z.xy*rmx;
	z.xyz=clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz;
	z*=2.0/clamp(dot(z.xyz,z.xyz),0.1,1.0);
	z+=p0;
	return length(z.xyz)/z.w;
}

float CircleOfConfusion(float t){//calculates the radius of the circle of confusion at length t
	return max(abs(focalDistance-t)*aperture,pixelSize*(1.0+t));
}
mat3 lookat(vec3 fw,vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}
float linstep(float a, float b, float t){return clamp((t-a)/(b-a),0.,1.);}// i got this from knighty and/or darkbeam
//random seed and generator
vec2 randv2;
float rand2(){// implementation derived from one found at: lumina.sourceforge.net/Tutorials/Noise.html
	randv2+=vec2(1.0,1.0);
	return fract(sin(dot(randv2 ,vec2(12.9898,78.233))) * 43758.5453);
}

float FuzzyShadow(vec3 ro, vec3 rd, float lightDist, float coneGrad, float rCoC){
	float t=0.01,d=1.0,s=1.0;
	for(int i=0;i<4;i++){
		if(t>lightDist)continue;
		float r=rCoC+t*coneGrad;//radius of cone
		d=DE(ro+rd*t)+r*0.66;
		s*=linstep(-r,r,d);
		t+=abs(d)*(0.8+0.2*rand2());
	}
	return clamp(s,0.0,1.0);
}

void main() {
	randv2=fract(cos((gl_FragCoord.xy+gl_FragCoord.yx*vec2(1000.0,1000.0))+vec2(time)*10.0)*10000.0);
	pixelSize=1.0/size.y;
	float tim=time*0.25;//camera, lighting and object setup
	float ct=cos(tim),st=sin(tim);
	rmx=mat2(ct,-st,st,ct);
	float z=cos(tim*0.3)*5.0;
	vec3 ro=vec3(vec2(ct,st)*(abs(z)+0.1)*(1.0+sin(tim*0.1)),z);
	vec3 rd=lookat(-ro,vec3(0.0,1.0,0.0)-sin(ro)*0.1)*normalize(vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,2.0));
  if (!setup_ray(eye, dir, ro, rd)) {  // boxplorify view
    return;
  }
	//focalDistance=min(length(ro)+0.001,1.0);
	//aperture=0.007*focalDistance;
	vec3 lightColor=vec3(1.0,0.5,0.25);
	vec4 col=vec4(0.0);//color accumulator, .w=alpha
	float t=0.0,mld=100.0;//distance traveled, minimum light distance
	for(int i=1;i<72;i++){//march loop
		if(col.w>0.9 || t>15.0)continue;//bail if we hit a surface or go out of bounds
		float rCoC=CircleOfConfusion(t);//calc the radius of CoC
		float d=DE(ro)+0.33*rCoC;
		float lightDist=DEL(ro);//the distance estimate to light
		mld=min(mld,lightDist);//the minimum light distance along the march
		if(d<rCoC){//if we are inside the sphere of confusion add its contribution
			vec3 p=ro-rd*abs(d-rCoC);//back up to border of SoC
			mcol=vec3(0.0);//clear the color trap
			bColoring=true;//collecting color samples with normal deltas
			vec2 v=vec2(rCoC*0.5,0.0);//use normal deltas based on CoC radius
			vec3 N=normalize(vec3(-DE(p-v.xyy)+DE(p+v.xyy),-DE(p-v.yxy)+DE(p+v.yxy),-DE(p-v.yyx)+DE(p+v.yyx)));
			bColoring=false;
			if(N!=N)N=-rd;//if no gradient assume facing us
			v=vec2(lightDist,0.0);//now find the light's general direction
			vec3 L=-normalize(vec3(-DEL(p-v.xyy)+DEL(p+v.xyy),-DEL(p-v.yxy)+DEL(p+v.yxy),-DEL(p-v.yyx)+DEL(p+v.yyx)));
			float lightStrength=1.0/(1.0+lightDist*lightDist*20.0);
			vec3 scol=mcol*0.1666*(0.2+0.4*(1.0+dot(N,L)))*lightStrength;//average material color * diffuse lighting * attenuation
			scol+=0.5*pow(max(0.0,dot(reflect(rd,N),L)),32.0)*lightColor*lightStrength;//specular lighting
			scol*=FuzzyShadow(p,L,lightDist,shadowCone,rCoC);//now stop the shadow march at light distance
			col.rgb+=lightColor/(0.5+mld*mld*5000.0)*(1.0-col.w);//add a bloom around the light
			mld=100.0;//clear the minimum light distance for the march
			float alpha=fudgeFactor*(1.0-col.w)*linstep(-rCoC,rCoC,-d);//calculate the mix like cloud density
			col=vec4(col.rgb+scol*alpha,clamp(col.w+alpha,0.0,1.0));//blend in the new color	
		}//move the minimum of the object and light distance
		d=abs(fudgeFactor*min(d,lightDist+0.33*rCoC)*(0.8+0.2*rand2()));//add in noise to reduce banding and create fuzz
		ro+=d*rd;//march
		t+=d;
	}//mix in background color
	vec3 scol=lightColor/(0.5+mld*mld*5000.0);//add one last light bloom
	col.rgb+=scol*(1.0-col.w);

	//gl_FragColor = vec4(clamp(col.rgb,0.0,1.0),1.0);
  write_pixel(dir, t, clamp(col.rgb, 0.0, 1.0));  // boxplorify write
}
