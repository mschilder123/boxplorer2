//Fur Trees by eiffie
//attempting some distance estimated fur and shading with one extra DE calc

#include "setup.inc"
#line 6

vec2 iResolution = vec2(xres, yres);
float iGlobalTime = time;

#define AUTO_OVERSTEP

//#define time iGlobalTime
#define size iResolution

#define TAO 6.283
vec2 rotate(vec2 v, float angle) {return cos(angle)*v+sin(angle)*vec2(v.y,-v.x);}
vec2 kaleido(vec2 v, float power){return rotate(v,floor(.5+atan(v.x,-v.y)*power/TAO)*TAO/power);}

vec2 kaleido6(vec2 v){return rotate(v,floor(0.5+atan(v.x,-v.y)*0.95493)*1.0472);}
vec2 kaleido12(vec2 v){return rotate(v,floor(0.5+atan(v.x,-v.y)*1.90986)*0.5236);}

vec3 mcol;//material color
mat2 r45=mat2(0.7071,0.7071,-0.7071,0.7071);
mat2 r30=mat2(0.866,0.5,-0.5,0.866);
mat2 rtrn=mat2(0.9689,-0.2474,0.2474,0.9689);

float DE(in vec3 z0){
	//mcol=vec3(0.5);
	//return length(z0)-1.0;
	if(z0.y>0.1)return z0.y+0.2;
	z0.xz=mod(z0.xz,2.0)-vec2(1.0);//+0.125*cos(floor(z0.zx*0.5)*2.0);
	float cyl=length(z0.xz);
	float d=100.0,dt=cyl+z0.y*0.025;
	for(int i=0;i<2;i++){
		vec3 z=z0;
		z.y-=float(i)*0.125;
		float c=floor(z.y*4.0);
		//z.yz=rotate(z.yz,-z.z*0.79*(1.0+c*0.1));
		float bm=-z.y-2.0+cyl*0.01;
		z.y=mod(z.y,0.25)-0.05;
		if(i==1)z.xz=z.xz*rtrn;
		z.xz=kaleido(z.xz,2.0-c);
		z.yz=rtrn*z.yz;
		bm=max(bm,-z.z+c*0.086);//0.065);
		dt=min(dt,max(max(abs(z.x),abs(z.y)),bm))-0.001-z.z*0.003;
		float c2=floor(z.z*16.0);
		z.z=mod(z.z,0.0625)-0.049;
		z.xy=rotate(z.xy,c2*0.25);
		z.xy=kaleido12(z.xy);
		z.yz=z.yz*r30;
		d=min(d,max(max(max(abs(z.x),abs(z.z)),-z.y-0.05+c*0.005),bm));
	}
	if(dt<d){
		d=dt;
		mcol=vec3(0.5,0.1,0.0);
	}else{
		mcol=vec3(0.5,0.6,0.2);
		mcol*=1.0+(-z0.y*0.75)*(z0.x+z0.z)/cyl;
	}
    mcol*=cyl + 0.5 + z0.y*0.5;//kind of what iq suggested
	return max(0.0,max(d,max(z0.y,-z0.y-2.0)));
}

float rndStart(vec2 co){return 0.1+0.9*fract(sin(dot(co,vec2(123.42,117.853)))*412.453);}

mat3 lookat(vec3 fw,vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,up));return mat3(rt,cross(rt,fw),fw);
}

void main(){
	float zoom=2.0,px=2.0/(size.y*zoom);//find the pixel size
	float tim=time*0.3;
	
	//position camera
	vec3 ro=vec3(0.5*sin(tim*0.43),-1.0,tim);
	vec3 rd=normalize(vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,zoom));
	rd=lookat(vec3(-1.15+0.75*sin(tim),-1.0,tim+2.0)-ro,vec3(0.0,1.0,0.0))*rd;
	//ro=eye;rd=normalize(dir);
	vec3 ld=normalize(vec3(0.4,0.75,0.4));//direction to light
	vec3 bcol=clamp(vec3(1.0-rd.y,1.0-rd.y-0.1*exp(-abs(rd.y*15.0)),1.0-0.5*exp(-abs(rd.y*5.0))),0.0,1.0);//backcolor
  if (!setup_ray(eye, dir, ro, rd)) return;
	//march
	
	float tG=abs((-2.0-ro.y)/rd.y),d,pd=10.0,os=0.0,step=0.0;
	vec2 g=ro.xz+rd.xz*tG;
	float t=DE(ro)*rndStart(gl_FragCoord.xy);
	float MIN_DIST=px*0.1;
	vec4 col=vec4(0.0);//color accumulator
	for(int i=0;i<78;i++){
		d=DE(ro+rd*t);
		float d1=max(d,px*t*0.5);
#ifdef AUTO_OVERSTEP
		if(d1>os){		//we have NOT stepped over anything
			if(t>tG)break;
			os=0.28*d1*d1/pd;//calc overstep based on ratio of this step to last
			step=d1+os;	//add in the overstep
			pd=d1;	//save this step length for next calc
		}else{
			step=-os;d1=1.0;pd=10.0;os=0.0;//remove ALL of overstep
		}
#else
			step=d1;
#endif
		if(d1<px*t){
			vec3 scol=mix(mcol,bcol,min(t*0.05,1.0));
			float d2=DE(ro+rd*t+ld*px*t);
			float shad=0.5*abs(d2/d);
			scol=scol*shad+vec3(0.2,0.05,-0.25)*(shad-0.5);
            vec3 pos = ro + t*rd;
            vec2 q = mod(pos.xz, 2.0) - vec2(1.0);
            scol *= clamp(length(q) + .8 + pos.y *.5,0.,1.);
			float alpha=(1.0-col.w)*clamp(1.0-d1/(px*t),0.0,1.0);
			col+=vec4(clamp(scol,0.0,1.0),1.0)*alpha;
			if(col.w>0.9)break;
		}
		t+=step;
	}

	//color the ground 
	if(rd.y<0.0){
		ro+=rd*tG;
		float s=1.0,dst=0.1;
		float t2=DE(ro)*rndStart(gl_FragCoord.xy);
		for(int i=0;i<4;i++){
			float d=max(0.0,DE(ro+ld*t2)*1.5)+0.05;
			s=min(s,3.0*d/t2);
			t2+=dst;dst*=2.0;
		}
		bcol*=0.8+0.2*s;
	}
	col.rgb+=bcol*(1.0-clamp(col.w,0.0,1.0));

	//gl_FragColor=vec4(col.rgb,1.0);
  write_pixel(dir, t, col);
} 
