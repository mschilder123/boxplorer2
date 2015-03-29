//Escaping Fleet by eiffie (another test of aa methods, not turned on by default)

#include "setup.inc"
#line 5

uniform int iters;
uniform int color_iters;
vec2 iResolution = vec2(xres, yres);

#define HI_QUAL
//#define USE_TEXTURE

#ifdef HI_QUAL
	#define FudgeFactor 0.75
	#define MarchSteps 96
	#define ShadowSteps 12
	#define AOSteps 5
#else
	#define FudgeFactor 1.0
	#define MarchSteps 64
	#define ShadowSteps 10
	#define AOSteps 4
#endif

//#define time iGlobalTime
#define size iResolution
#define tex iChannel0
const vec3 LightDir=vec3(0.0,1.0,0.0),LightColor=vec3(1.0,0.99,0.9),SkyColor=vec3(0.1,0.16,0.27);

float PixelSize,st;
vec3 mcolor;

float rnd(vec2 co){return fract(sin(dot(co,vec2(123.42,117.853)))*412.453);}
float noyz(vec2 p){//iq's hash and noise functions
	vec2 c=floor(p),f=fract(p),v=vec2(1.0,0.0);
	return mix(mix(rnd(c),rnd(c+v.xy),f.x),mix(rnd(c+v.yx),rnd(c+v.xx),f.x),f.y);
}
vec3 texture(vec2 p){
#ifdef USE_TEXTURE
	return texture2D(tex,p*10.).rgb*2.0;
#else
	float n=noyz(p*256.0);
	return vec3(0.8)+vec3(0.4*n,0.0,0.7-n);
#endif
}

const float mr=0.13, SCALE = -1.7;
vec4 scale=vec4(SCALE,SCALE,SCALE,abs(SCALE));
const vec4 p0=vec4(2.0,-0.32,2.48,1.0);
float DE(in vec3 z0){//amazing box by tglad
	vec4 z = vec4(z0,1.0),zG;
	for (int n = 0; n < iters; n++) {
		z.xyz=clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz;
		z*=scale/clamp(dot(z.xyz,z.xyz),mr,1.0);
		if(n==3)zG=z;
		z+=p0;
	}
	float dG=(length(max(abs(zG.xyz)-vec3(0.8,4.2,0.0),0.0))-0.01)/zG.w;
	return min(dG,(length(max(abs(z.xyz)-vec3(4.4,0.9,1.5),0.0))-0.01)/z.w);
}
float CE(in vec3 z0){//same for coloring
	vec4 z = vec4(z0,1.0),zG;
	for (int n = 0; n < color_iters; n++) {
		z.xyz=clamp(z.xyz, -1.0, 1.0) *2.0-z.xyz;
		z*=scale/clamp(dot(z.xyz,z.xyz),mr,1.0);
		if(n==3)zG=z;
		z+=p0;
	}
	float dG=length(max(abs(zG.xyz)-vec3(0.8,4.2,0.0),0.0))/zG.w;
	float dS=length(max(abs(z.xyz)-vec3(4.4,0.9,1.5),0.0))/z.w;
	vec3 col=texture(10.0*z0.xy+5.0*z0.zz);
	if(dS<dG){dS-=col.r*0.01/z.w;}
	else {col=col.brg;dG-=col.r*0.01/zG.w;col.r+=0.5;}
	mcolor+=col;
	return min(dS,dG);
}

//I remember someone doing a similar background, sorry I forgot who.
vec3 getBackground( in vec3 rd ){
	vec3 bcol=SkyColor+rd*0.1+LightColor*max(0.0,dot(rd,LightDir))*0.2;
	float y=1.0-abs(rd.y),a=0.44+atan(rd.x,rd.z);
	vec2 pt=vec2(a+sin(7.0*y+a*10.0+time*0.25)*0.05*y,rd.y+time*0.25);
	bcol*=texture(pt);
	bcol+=LightColor*pow(max(0.0,dot(rd,LightDir)),110.0+st*100.0);
	return bcol;
}

float shadao(vec3 ro, vec3 rd, float px, float max_dist){//pretty much IQ's SoftShadow
	float res=1.0,d,t=2.0*px,min_step=0.5*PixelSize;
	for(int i=0;i<ShadowSteps;i++){
		d=max(0.0,DE(ro+rd*t))+min_step;
		if(t+d>max_dist)break;
		t+=d;
		res=min(res,3.0*d/t);
	}
	return res;
}

float fakeAO(vec3 ray, vec3 norm, float ao_eps) {//from rrrola
	float ao=1.0,w=0.1/ao_eps,dist=2.0*ao_eps,d;
	for (int i=0; i<AOSteps; i++) {
		d = DE(ray + norm*dist);
		ao -= (dist-d) * w; 
		w *= 0.5; dist = dist*2.0 - ao_eps;
	}
  	return clamp(ao, 0.0, 1.0);
}

vec3 shade(in vec3 ro, in vec3 rd, in float t, in vec3 color){
	float px=PixelSize*t;
	vec2 v=vec2(0.5*px,0.0);
	float ds=DE(ro+rd*t);
	ro+=rd*(t+ds-px);
	mcolor=vec3(0.0);//clear material color before taking samples
	float d=CE(ro);
	vec3 dn=vec3(CE(ro-v.xyy),CE(ro-v.yxy),CE(ro-v.yyx));
	vec3 dp=vec3(CE(ro+v.xyy),CE(ro+v.yxy),CE(ro+v.yyx));
	mcolor*=0.143;
	vec3 norm=normalize(dp-dn);
	float shad=shadao(ro,LightDir,px,10.0);
	float ao=fakeAO(ro,norm,px);
	float dif=dot(norm,LightDir)*0.5+0.5;
	vec3 reflDir=reflect(rd,norm);
	float refl=dot(LightDir,reflDir);
	dp=abs(vec3(d)-0.5*(dp+dn));
	d=max(dp.x,max(dp.y,dp.z))/(0.04*px);//calc curvature to remove sparkles
	//float spec=mspec*pow(max(0.0,refl),mspecExp)*clamp(1.0-d,0.0,1.0)*mcolor.r;
    float spec=max(0.0,refl)*clamp(1.0-d,0.0,1.0)*mcolor.r;
	dif=min(dif,shad+0.1);
	vec3 diffuse_col=mcolor+vec3(0.1,0.0,-0.1)*refl;
#ifdef HI_QUAL
	vec3 bcol=getBackground(reflDir);
#else
	vec3 bcol=LightColor*pow(max(0.0,dot(reflDir,LightDir)),110.0+st*100.0);
#endif
	return mix(ao*(diffuse_col*dif+bcol*shad*spec)*(1.25-0.25*st),color,clamp(ds/px,0.0,1.0));
}

vec3 scene( vec3 ro, vec3 rd, out float total )
{//march
	float t=DE(ro)*rnd(gl_FragCoord.xy)*0.75;
	float d,dm=100.0,tm=0.0,MIN_DIST=PixelSize*0.001,od=1000.0;
	bool bGrab=false;
	vec4 hit=vec4(-1.0);
	for(int i=0;i<MarchSteps;i++){
		d=DE(ro+rd*t)*FudgeFactor;
 #ifdef HI_QUAL
		if(d>od){
			if(bGrab && od<PixelSize*(t-od) && hit.x<0.0){
				hit.x=t-od;
				hit=hit.yzwx;
				bGrab=false;
			}
		}else bGrab=true;
		od=d;
 #endif
		if(d<dm){tm=t;dm=d;}//save the max occluder
		t+=d;
		if(t>10.0 || d<MIN_DIST)break;
	}
  total = t;
	//we have saved the edges but there is also a min or final surface
	if(tm>hit.w && dm<PixelSize*(tm-dm)){//if minimum has not been saved
		if(hit.x>0.0)hit=hit.wxyz;//write over the last entry, not the first
		hit.x=tm-dm;hit=hit.yzwx;
	}else if(od<PixelSize*(t-d) && hit.x<0.0){//save final distance if we have room on the stack
		hit.x=t-d;hit=hit.yzwx;	
	}

	//color the background 
	vec3 col=getBackground(rd);
	
	//add in the object(s)
#ifdef HI_QUAL
	for(int i=0;i<4;i++){//play back the hits and mix the color samples
#endif
		hit=hit.wxyz;//pop (when hi qual is off we are just using the max occluder)
		if(hit.x>0.0)col=shade(ro,rd,hit.x,col);
#ifdef HI_QUAL
	}
#endif
	
	return clamp(col,0.0,1.0);
}
	 
mat3 lookat(vec3 fw,vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}

void main() {
	st=sin(time*2.5);
	vec3 ro=vec3(sin(time*0.02)*2.31,-0.005,0.051);//I should have put more thought into this :)
	vec3 dr=vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,2.0);
	vec3 rd=normalize(dr);
	PixelSize=2.5/(size.y*dot(rd,dr));
	vec3 fw=mix(vec3(sin(time*0.1),cos(time*0.27),sin(time*0.13)),-ro,clamp(dot(ro,ro)*0.12,0.0,1.0));
	float d=DE(ro);
	fw=mix(vec3(sign(sin(time*0.02+1.57)),0.0,0.0),fw,smoothstep(-0.025,0.04,d*d));
	rd=lookat(fw,vec3(0.3+0.5*sin(time*0.23),1.0,0.4))*rd;
  if (!setup_ray(eye, dir, ro, rd)) return;
  float total;
	vec3 color=scene(ro,rd,total);
	if(color!=color)color=vec3(0.0,1.0,0.0);
  write_pixel(dir, total, color);
} 
