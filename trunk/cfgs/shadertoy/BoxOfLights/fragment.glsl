// Box-O-Lights by eiffie
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#include "setup.inc"
#line 5
vec2 iResolution = vec2(xres, yres);
uniform int iters;

#define LIGHT_AA
//#define HIQUAL
const int MarchSteps=32,ShadowSteps=24;//HIQUAL is a bit misleading with 32 march steps!

//#define time iGlobalTime
#define size iResolution
float HitDistance=0.666/size.y;
const float shadows=0.75,contrast=0.75,spec=1.0,specExp=32.0,maxDepth=10.0;
const vec3 lightColor=vec3(1.5,0.99,0.9);
bool bColoring=false,bCheckLight=true;
vec3 dc;//some wonderful globals with opaque names :)
float mld;
float L1=cos(time*2.0)*0.3,L2=sin(time*2.0);
vec3 posL=vec3(0.5+L1,0.5-L1,L2);//position of lights

const float mr=-0.26, mxr=0.74;
const vec4 scale=vec4(-2.36,-2.36,-2.36,2.36);
const vec3 p0=vec3(0.72,3.44,-1.6),pS=vec3(0.25,0.25,0.15),rc=vec3(2.1,10.0,10.0);
float DEL(in vec3 z){//distance estimate to light
	z=clamp(z, -1.0, 1.0) *2.0-z;
	z.xy=abs(z.xy);
	z*=scale.x;
	return (length(z+posL))/scale.w;
}
float DE(in vec3 z0){//amazing box by tglad with mods
	vec4 z = vec4(z0,1.0),p2=vec4(p0+abs(z0)*pS,0.15);
	float dL=100.0;
	if(bCheckLight){dL=DEL(z0);mld=min(mld,dL);}
	for (int n = 0; n < iters; n++) {
		z.xyz=abs(clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz);
		if(z.x<z.y)z.xy=z.yx;
		z*=scale/clamp(dot(z.xyz,z.xyz)+0.046,mr,mxr);
		z+=p2;
	}
	//float dS=min(dL,(max(abs(z.x),max(abs(z.y),abs(z.z)))-ABSCL)/z.w);
	float dS=min(dL,(length(max(abs(z.xyz)-rc,0.0))-1.0)/z.w);
	if(bColoring){
		if(dS==dL)dc=vec3(1.0);
		dc+=z.xyz;
	}
	return dS;
}

float shadow(in vec3 ro, in vec3 rd, float max_depth, float k)
{
	float fStep=1.0,t=HitDistance*4.0,d=1.0;	
	for(int i=0;i<ShadowSteps;i++){
		if(t>max_depth || d<HitDistance)break;
		t+=d=DE(ro+rd*t);
		fStep+=1.0;
	}
	if(d<HitDistance)return 1.0;
	return min((fStep-d*k/t)/float(ShadowSteps),1.0);
}

vec3 scene( vec3 ro, vec3 rd, out float total )
{// find color of scene
	vec3 col=vec3(0.0);
	float t=0.0,d=1.0,fStep=1.0;
	mld=100.0;bCheckLight=true;
	for(int i=0;i<MarchSteps;i++){
		if(t>maxDepth || d<HitDistance)break;
		t+=d=DE(ro+rd*t);
		fStep+=1.0;
	}
	if( d< 0.25 ){//hit something	
		if(dc.b==1.0){//hit a light
			col=lightColor;
		}else{
        		ro+= rd * (t-HitDistance);// advance ray position
			bColoring=true;bCheckLight=false;
			vec2 ve=vec2(HitDistance,0.0);
			dc=vec3(0.0);
			float d=DE(ro),d1=DE(ro-ve.xyy),d2=DE(ro+ve.xyy);
			float d3=DE(ro-ve.yxy),d4=DE(ro+ve.yxy);
			float d5=DE(ro-ve.yyx),d6=DE(ro+ve.yyx);
			bColoring=false;
			vec3 diffuse=vec3(0.5)+sin(dc*0.143)*0.4;
			ve=vec2(DEL(ro),0.0);
			vec3 lightDir=-normalize(vec3(-DEL(ro-ve.xyy)+DEL(ro+ve.xyy),
				-DEL(ro-ve.yxy)+DEL(ro+ve.yxy),-DEL(ro-ve.yyx)+DEL(ro+ve.yyx)));
			float shad=shadow( ro, lightDir, ve.x, 1.0 );
			float lightStrength=(1.0-shadows*shad)/max(ve.x*ve.x*12.0,0.5);
#ifndef LIGHT_AA
			vec3 nor=normalize(vec3(-d1+d2,-d3+d4,-d5+d6));
			vec2 lacc=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
#else //we are despeckling the specular highlights :)
			vec3 nor=normalize(vec3(-d1+d,-d3+d,-d5+d));if(nor!=nor)nor=-rd;
			vec2 lacc=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			nor=normalize(vec3(-d+d2,-d+d4,-d+d6));if(nor!=nor)nor=-rd;
			lacc+=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			nor=normalize(vec3(-d+d2,-d3+d,-d5+d));if(nor!=nor)nor=-rd;
			lacc+=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			nor=normalize(vec3(-d+d2,-d+d4,-d5+d));if(nor!=nor)nor=-rd;
			lacc+=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			nor=normalize(vec3(-d1+d,-d+d4,-d5+d));if(nor!=nor)nor=-rd;
			lacc+=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			nor=normalize(vec3(-d1+d,-d3+d,-d+d6));if(nor!=nor)nor=-rd;
			lacc+=vec2(max(0.0, dot(lightDir, nor)),pow(max(0.0,dot(rd,reflect(lightDir,nor))),specExp));
			lacc*=0.1666666;//lighting calculations averaged across each facet
#endif
			col=lightStrength*(diffuse*(1.0-contrast*(1.0-lacc.x))+lightColor*spec*lacc.y);
			
		}
	}
  total = t;
	col=mix(col,abs(rd.xzy)*0.05,pow(fStep/float(MarchSteps),4.0));
	col=max(col,lightColor/max(mld*mld*10000.0,1.0));
	return vec3(clamp(col,0.0,1.0));
}	 

mat3 lookat(vec3 fw,vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}

void main() {
	vec2 v=vec2(cos(time*0.2),sin(time*0.3)*1.5);
	vec3 ro=vec3(v*(0.25+v.x),v.y*0.47-0.97);
	mat3 rotCam=lookat(-ro,vec3(-1.0,0.0,0.0));
	vec3 rd=rotCam*normalize(vec3((2.0*(gl_FragCoord.xy)-size.xy)/size.y,1.0));
  float total;
  if (!setup_ray(eye, dir, ro, rd)) return;
	vec3 col=scene(ro,rd, total);
#ifdef HIQUAL
	rd=rotCam*normalize(vec3((2.0*(gl_FragCoord.xy+vec2(0.5,0.0))-size.xy)/size.y,1.0));
	col+=scene(ro,rd);
	rd=rotCam*normalize(vec3((2.0*(gl_FragCoord.xy+vec2(0.0,0.5))-size.xy)/size.y,1.0));
	col+=scene(ro,rd);
	rd=rotCam*normalize(vec3((2.0*(gl_FragCoord.xy+vec2(0.5))-size.xy)/size.y,1.0));
	col+=scene(ro,rd);
	col*=0.25;
#endif
//	gl_FragColor = vec4(col,1.0);
  write_pixel(dir, total, col);
}
