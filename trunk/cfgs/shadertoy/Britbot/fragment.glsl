//Britbot by eiffie
//A little bit like mode7 from poljere https://www.shadertoy.com/view/ltsGWn
//...and a lot of tunnel from iq https://www.shadertoy.com/view/Ms2SWW

#include "setup.inc"
#line 7

#define iGlobalTime time
vec2 iResolution = vec2(xres, yres);
#define iChannel0 iBackbuffer

#define PI 3.14159
vec3 tunnel(vec3 rd){
	vec2 uv=rd.xy/rd.z;
	float tm=iGlobalTime*sign(rd.z)*1.5;
	float pw=0.5+pow(min(abs(sin(tm*0.1))+0.25,1.0),16.0)*16.0;
	float r=pow(pow(uv.x*uv.x,pw)+pow(uv.y*uv.y,pw),0.5/pw);
	float x;
	for(int i=1;i<10;i++){//what kind of maths are these?
		x=uv.x+sin((0.5/r+0.5*tm)*2.0)*float(i)*float(i)*0.001;
        pw=0.5+pow(min(abs(sin((tm+max(1.0-r,0.0))*0.1))+0.25,1.0),16.0)*16.0;
		r=pow(pow(x*x,pw)+pow(uv.y*uv.y,pw),0.5/pw);
	}
	uv.x=x;
	float a=atan(uv.y,uv.x)/3.14159;
	vec2 p=vec2(0.5/r+0.5*tm,a)*8.0;
	p.y*=sign(uv.x);
	vec2 c=floor(p);
	p=fract(p);
	p.x=pow(p.x,clamp(abs(rd.z)+pw/16.0,1.0,2.0));
	uv=p;
	uv=2.0*(uv-0.5);
	float r2=pow(pow(uv.x*uv.x,pw)+pow(uv.y*uv.y,pw),0.5/pw);
	p=clamp(p*1.5-0.25,0.0,1.0);
	vec3 col=vec3(0.5)+0.5*sin(vec3(c.xy,c.x+c.y));
    if((a<0.25 && a>-0.25) || a<-0.75 || a>0.75){
        float d=max(abs(p.x-0.5),abs(p.y-0.5))-0.5;
        if(d<0.0){//min(p.x,p.y)>0.0 && max(p.x,p.y)<1.0){
        	if(rd.z<0.0)p.y=1.0-p.y;
			col=mix(col,texture2D(iChannel0,p).rgb,smoothstep(0.0,0.05,-d));
        }
    }
	col*=2.0*pow(r,1.75)*clamp(3.0-r2*3.0,0.0,1.0);
	if(col!=col)col=vec3(0.0);
	return clamp(col,0.0,1.0);
}
mat3 lookat(vec3 fw, vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}
struct intersect{float t, d; vec3 N;}I1,I2; 
void zStack(intersect I, float px){
	if(I.t<=0.0 || I.d>px*I.t)return;
	if(I.t<I1.t){I2=I1;I1=I;}
	else if(I.t<I2.t)I2=I;
}
#define maxDepth 10.0
//pS=p1-ro, pD=p2-p1 Hopefully no one thinks this actually works!
intersect Segment(in vec3 pS, in vec3 pD, in float r, in vec3 rd){//mod from iq's
	intersect intr=intersect(0.0,maxDepth,vec3(0.0));
	float d=dot(rd,pD);
	float t=clamp((dot(rd,pS)*d-dot(pS,pD))/(dot(pD,pD)-d*d),0.0,1.0);
	pS+=pD*t;
	intr.N=-pS;
	float b=dot(pS,rd);
	float h=b*b-dot(pS,pS);
	d=sqrt(max(0.0,-h))-r;
	intr.d=max(0.0,d);
	intr.t=b+min(d,0.0)-sqrt(max(0.0,h+r*r));
	return intr;
	//dist: intr.t
	//aa/dof: clamp(intr.d/(px*intr.t),0.0,1.0);
	//shad: clamp(k*intr.d/intr.t,0.0,1.0);
	//normal: normalize(rd*intr.t+intr.N);
}
vec3 Light(intersect I, vec3 rd, float px, vec3 col){
	float aac=1.0-clamp(I.d/(px*I.t),0.0,1.0);
	if(aac>0.0){
		vec3 N=normalize(rd*I.t+I.N);
		vec3 L=normalize(vec3(0.5,0.8,0.4));
		vec3 R=reflect(rd,N);
		col=mix(col,(vec3(1.0,0.3,0.4)+0.2*tunnel(R))*(0.5+0.5*dot(L,N)),aac);
	}
	return col;
}
vec3 jsolve( vec3 a, vec3 b, float ln, vec3 rt )//mod from iq's
{//simple joint with equal lengths
	vec3 p=b-a,q=p*0.5;
	return a+q+sqrt(max(0.0,ln*ln-dot(q,q)))*normalize(cross(p,rt));
}
vec3 britbot(vec3 ro, vec3 rd, vec3 col){
	float px=2.5/iResolution.y,tm=iGlobalTime*10.0;
	I1.t=I2.t=I1.d=I2.d=maxDepth;
	float ct=cos(tm),st=sin(tm),st2=sin(tm*0.3);
	float h=(ct+st)*ct*-0.25;
	vec3 b1=vec3(0.0,h,0.0),b2=vec3(st2*0.2,-0.75+h,-0.2);
	
	vec3 le,lh=vec3(-0.75+0.1*st,-0.4-0.2*ct,0.3-0.3*ct),re,rh=vec3(0.75+0.1*st,-0.4-0.2*st,0.3+0.3*ct);
	vec3 lk,lf=vec3(-0.25,-2.0+max(0.0,ct*0.5),-0.25+0.4*st),rk,rf=vec3(0.25,-2.0+max(0.0,-ct*0.5),-0.25-0.4*st);
	vec3 rt=normalize(vec3(1.0+0.4*st2,-0.4*st2,0.0));
	le=jsolve(b1,lh,0.6,rt.yxz);
	re=jsolve(b1,rh,0.6,-rt.yxz);
	lk=jsolve(b2,lf,0.7,rt);
	rk=jsolve(b2,rf,0.7,rt);
	float TR=0.15;
	vec3 n=vec3(TR,0.0,-TR*0.5);
	zStack(Segment(b1-ro,b2-b1,TR,rd),px);
	TR*=0.75;
	zStack(Segment(lk-ro,b2-n.xyy-lk,TR,rd),px);
	zStack(Segment(rk-ro,b2+n.xyy-rk,TR,rd),px);
	TR*=0.75;
	zStack(Segment(lf-ro,lk-lf,TR,rd),px);
	zStack(Segment(rf-ro,rk-rf,TR,rd),px);
	TR*=0.75;
	zStack(Segment(le-ro,b1-n.xyy-le,TR,rd),px);
	zStack(Segment(re-ro,b1+n.xyy-re,TR,rd),px);
	TR*=0.75;
	zStack(Segment(lh-ro,le-lh,TR,rd),px);
	zStack(Segment(rh-ro,re-rh,TR,rd),px);
	TR*=0.75;
	mat3 mx=lookat(vec3(st2*ct,0.25*(st2+st),1.0),vec3(0.25*st2*st,1.0,0.25*st2*ct));
	vec3 n1=vec3(0.0,0.25+h,0.0)*mx,h1=vec3(-0.4,0.25+h,0.0)*mx,h2=vec3(-0.4,0.75+h,0.0)*mx;
	vec3 h3=vec3(0.4,0.75+h,0.0)*mx,h4=vec3(0.4,0.25+h,0.0)*mx;
	
	n=vec3(0.0,0.0,1.0)*mx;
	float t=-dot(n,ro)/dot(n,rd);
	if(t>0.0){
		vec3 p=mx*(ro+rd*t);
		p.y-=h;
		if(p.x>-0.4 && p.x<0.4 && p.y>0.25 && p.y<0.75){
			col=texture2D(iChannel0,vec2(1.25,2.0)*(p.xy+vec2(0.4,-0.25))).rgb;
		}
	}
	zStack(Segment(b1-ro,n1-b1,TR,rd),px);
	zStack(Segment(h1-ro,h2-h1,TR,rd),px);
	zStack(Segment(h2-ro,h3-h2,TR,rd),px);
	zStack(Segment(h3-ro,h4-h3,TR,rd),px);
	zStack(Segment(h1-ro,h4-h1,TR,rd),px);
	
	col=Light(I2,rd,px,col);
	col=Light(I1,rd,px,col);
	return col;
}

void main(void)
{
	vec2 uv = (2.0*gl_FragCoord.xy-iResolution.xy)/ iResolution.y;
	float tm=iGlobalTime*0.1;
	vec3 ro=vec3(sin(tm),0.0,cos(tm)*3.0);
	tm=abs(sin(tm*1.5));
	vec3 up=vec3(1.0-tm,1.0,0.0);
	vec3 rd=lookat(-ro,up)*normalize(vec3(uv,1.0));
  if (!setup_ray(eye, dir, ro, rd)) return;
 	vec3 col=tunnel(rd);
	col=britbot(ro,rd,col);
	gl_FragColor = vec4(col,1.0);
}
