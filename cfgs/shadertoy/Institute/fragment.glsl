//The Institute by eiffie

#include "setup.inc"
#line 5

vec2 size = vec2(xres, yres);

#define AUTO_OVERSTEP

float px;
vec4 prp=vec4(0.0);
vec3 L;

float rnd(vec2 c){return fract(sin(dot(vec2(1.317,19.753),c))*413.7972);}
float rndStart(){
	return 0.5+0.5*rnd(gl_FragCoord.xy);
}
float smin(float a,float b,float k){float h=clamp(0.5+0.5*(b-a)/k,0.0,1.0);return b+h*(a-b-k+k*h);}
float noyz(vec2 p){//value noise - original author??
	vec2 c=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);
	float k=257.0,u=c.x+c.y*k;
	vec4 v=vec4(u,u+1.0,u+k,u+k+1.0);v=fract(fract(v*1.23456789)*9.18273645*v);
	return mix(mix(v.x,v.y,f.x),mix(v.z,v.w,f.x),f.y);
}
float fbm(vec2 p){
	vec2 s=sin(p*0.3+2.4*sin(p.yx*0.3));
	float h=1.0+(s.x+s.y)*0.5,a=0.5;
	p+=s.yx*0.66;
	for(int i=0;i<3;i++){
		h+=noyz(p)*a;
		a*=0.3;p=vec2(p.x+p.y*0.7,p.x-p.y+0.13)*2.0;
	}	
	return h;
}
vec3 Sky(vec3 rd){
	vec3 col=mix(vec3(0.4,0.2,0.0),vec3(0.5,0.7,0.7),clamp(0.25-0.5*rd.z,0.0,1.0));
	col+=vec3(1.2,1.1,0.9)*pow(max(0.0,dot(rd,L)),100.0);
	float h=noyz(rd.xy*5.0)*0.5+noyz(vec2(2.0*abs(atan(rd.y,rd.x)),rd.z*20.0));
	col=mix(col,vec3(0.9,0.8,0.6),clamp(h*0.5-0.4,0.0,1.0));
	return mix(vec3(0.1,0.2,0.1),col,clamp((-rd.z+0.05)*20.0,0.0,1.0));
}

float capsule(vec3 p){return length(vec3(p.x,p.y-clamp(p.y,-3.0,10.0),p.z));}
float rcap(vec3 p){return length(vec2(p.x-clamp(p.x,-7.0,7.0),max(abs(p.y),abs(p.z))));}

float DE(in vec3 p){
	float outlay=(length(p.xy+vec2(8.0,-6.0))-45.0);
	outlay=smin(outlay,abs(p.y-28.0+sin(p.x*0.01)*40.0)-5.0,20.0)*0.02;
	float h=fbm(p.xy*0.1)*10.0*pow(clamp(outlay,0.0,1.0),0.7),d=-p.z+7.0-h,dG=d;
	for (int n = 0; n < 4; n++) {
		p=clamp(p, -3.1, 3.1) *2.0-p;
		p+=vec3(1.9,3.2,8.6);
		p=p.yzx;
		d=min(d,min(capsule(p)-3.0,rcap(p)-2.0));
	}
	
	if(d<0.25 && d<dG){
		float flr=floor(p.y);
		float rs=1.0+sin(flr)*0.5;
		p=vec3(rs,1.0,rs)*0.5-abs(mod(p,vec3(rs,1.0,rs))-vec3(rs,1.0,rs)*0.5);
		float d2=d+0.05;
		if(flr<-rs || rs<0.6)
			d=min(d2,max(d,min(p.y,min(p.x,p.z))-0.0125));
		else if(flr<rs*2.0 || rs<0.75)
			d=max(d,min(p.y,min(p.x,p.z))-0.025);
		else {
			d=max(d,min(p.y,max(p.x,p.z))-0.05);
			if(prp.x<0.0)prp=vec4((rs>1.3)?length(p):10.0,0.0,0.0,(max(p.x,p.z)<0.06)?4.0:3.0);
		}
		if(prp.x<0.0){
			if(d==d2)prp=vec4(10.0,0.0,0.0,2.0);
			else prp=vec4(10.0,0.0,0.0,3.0);
		}
	}else if(prp.x<0.0)prp=vec4(10.0,h,outlay,1.0);
	return d;
}

float shadao(vec3 ro, vec3 rd, float px){//pretty much IQ's SoftShadow
	float res=1.0,d,t=10.0*px*rndStart();
	for(int i=0;i<12;i++){
		d=max(0.0,DE(ro+rd*t)*1.5);
		t+=d;
		res=min(res,d/t);
	}
	return res;
}

vec3 Color(vec3 ro, vec3 rd, float t, vec3 col, bool bFill){
	ro+=rd*t;
	prp.x=-1.0;
	float d=DE(ro),spec=0.0,n=noyz(ro.xy*3.0);
	vec2 e=vec2(px*t,0.0);
	vec3 dn=vec3(DE(ro-e.xyy),DE(ro-e.yxy),DE(ro-e.yyx));
	vec3 dp=vec3(DE(ro+e.xyy),DE(ro+e.yxy),DE(ro+e.yyx));
	vec3 N=(dp-dn)/(length(dp-vec3(d))+length(vec3(d)-dn));
	vec3 R=reflect(rd,N);
	vec3 lc=vec3(1.0,0.9,0.8),sc,rc=Sky(R);
	if(prp.w<1.5){
		sc=mix(vec3(0.8,0.8,0.3),vec3(0.2),clamp(abs(abs(ro.y-28.0+sin(ro.x*0.01)*40.0)-0.25)*4.0,0.0,1.0));
		sc=mix(sc,mix(vec3(0.2,0.3,0.1),vec3(0.6,0.5,0.4),clamp(prp.y*0.05-0.25,0.0,1.0)),clamp(prp.z*100.0,0.0,1.0));
		spec=0.2;//clamp(2.5-prp.y,0.0,1.0);
		n*=0.1;
	}else if(prp.w<2.5){sc=vec3(0.4,0.5,0.6);spec=1.0;n*=0.3;
	}else if(prp.w<3.5){sc=vec3(0.6,0.63,0.62);spec=0.5;n*=0.4;}
	else {sc=vec3(0.5,0.5,0.0);spec=0.5;}
	sc*=(1.0-n);
	float sh=clamp(shadao(ro,L,px*t)+0.2,0.0,1.0);
	sh=sh*(0.5+0.5*dot(N,L))+exp(-prp.x*5.0);
	vec3 scol=sh*lc*(sc+0.5*spec*rc*pow(max(0.0,dot(R,L)),2.0));
	if(bFill)d*=0.05;
	col=mix(scol,col,clamp(d/(px*t),0.0,1.0));
	return col;
}
mat3 lookat(vec3 fw){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,vec3(0.0,0.0,-1.0)));return mat3(rt,cross(rt,fw),fw);
}
vec3 cpnt(float t){
	if(t<0.5)return vec3(-45.0,30.0,-38.0);
	if(t<1.5)return vec3(-20.0,23.0,-25.0);
	if(t<2.5)return vec3(-13.0,23.0,-2.5);
	if(t<3.5)return vec3(10.0,26.0,-2.5);
	return vec3(12.0,21.0,6.6);
}
vec3 path(float t){
	float t2=t-t*t*0.05;
	if(t<10.0)return vec3(-1000.0+t2*200.0,30.0-sin((-1000.0+t2*200.0)*0.01)*40.0,4.25);
	t2=time-10.0;
	float r=60.0-t2;
	if(t<25.0)return vec3(-8.0+r*cos(t2*0.3),6.0+r*sin(t2*0.25),-t2*4.0);
	if(t<45.0){
		t2=(t-25.0)/5.0;r=floor(t2);t2=fract(t2);
		return mix(cpnt(r),cpnt(r+1.0)-0.25*cpnt(r+2.0)*(1.0-t2),t2);
	}
	if(t<60.0){
		t2=t-45.0;
		t2=t2-t2*t2*0.025;
		r=500.0-t2*50.0;
		vec3 p=vec3(cos(t2*0.2)*r,sin(t2*0.2)*r,-30.0);
		float d1=DE(p),d2=DE(p+vec3(0.0,-10.0,0.0));
		p.z+=(d1+d2)*0.4;
		return p;
	}
	return vec3(-8.0+cos(t*0.3)*60.0,6.0+60.0*sin(t*0.2),-29.5+sin(t*0.1)*10.0);
}
void main() {
	px=0.5/size.y;
	L=normalize(vec3(0.5,0.3,-0.6));
#if 0
	float tim=time;
	vec3 ro=path(tim);
	vec3 ta=path(tim+0.5);
	if(tim>10.0 && tim<25.0)ta=vec3(0.0);
	if(tim>60.0)ta=vec3(0.0,0.0,-20.0);
	vec3 rd=lookat(ta-ro)*normalize(vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,3.0));
#else
  vec3 ro, rd;
#endif
  if (!setup_ray(eye, dir, ro, rd)) return;
	//ro=eye*10.0;rd=normalize(dir);
	float t=DE(ro)*rndStart(),d=0.0,od=1.0,step=0.0,os=0.0,pd=10.0;
	if(t<0.0){
		gl_FragColor=vec4(0.0,0.0,0.0,1.0);
		return;
	}
	vec4 edge=vec4(-1.0);
	bool bGrab=false;
	for(int i=0;i<78;i++){
		t+=step;
		d=DE(ro+rd*t);
#ifdef AUTO_OVERSTEP
		if(d>=os){		//we have NOT stepped over anything
			os=0.36*d*d/pd;//overstep based on ratio of this step to last
			step=d+os;	//add in the overstep
			pd=d;		//save this step length for next calc
		}else{step=-os;d=0.0001;pd=100000.0;os=0.0;}//remove overstep
#else
		step=d;
#endif

		if(d>od){
			if(bGrab && od<px*t && edge.x<0.0){
				edge=vec4(edge.yzw,t-od);
				bGrab=false;
			}
		}else bGrab=true;
		od=d;
		if(t>1000.0 || d<0.00001)break;
	}
	bool bFill=false;
	d*=0.05;
	if(d<px*t){
		if(edge.x>0.0)edge=edge.wxyz;
		edge=vec4(edge.yzw,t);
		bFill=true;
	}
	vec3 col=Sky(rd);
	for(int i=0;i<4;i++){
		if(edge.w>0.0)col=Color(ro,rd,edge.w,col,bFill);
		edge=edge.wxyz;
		bFill=false;
	}
	//float dimmer=clamp(min(abs(time-10.0),min(abs(time-25.0),min(abs(time-45.0),abs(time-60.0)))),0.0,1.0);
	//gl_FragColor = vec4(1.5*col*dimmer*dimmer,1.0);
  write_pixel(dir, t, 1.5*col);
}
