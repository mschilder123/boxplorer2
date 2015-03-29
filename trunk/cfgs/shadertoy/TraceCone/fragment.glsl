//https://www.shadertoy.com/view/XsBXRV
//precision highp float;

#include "setup.inc"
#line 6

#define iGlobalTime time
vec2 iResolution = vec2(xres, yres);
vec4 iMouse = vec4(0.);

#define float3 vec3
#define float2 vec2
#define float4 vec4
#define float3x3 mat3

float3 campos=float3(-10.0,2.0,0.0);
float3 look_at=float3(0.0,1.0,0.0);
float3 up=float3(0,1,0);
float3 forward;
float3 right;

float3 light=float3(0,10,10);

const float MAX_RAY_LENGTH=10000.0;

float3 pixelRay(float2 uv)
{
	uv*=0.75;
    return normalize(forward+up*uv.y+right*uv.x);
}

void RP(float3 tp0, float3 dp1, float3 dp2, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
	float3 dp0=rp0-tp0;

	float3 dett =cross(dp1,dp2);
	float3 detuv=cross(dp0,rd);

	float det=(-1.0)/dot(dett,rd);

	float u=(dot(detuv,dp2))*det;
	float v=(dot(detuv,dp1))*det;
	t=(dot(dett ,dp0))*det;
    if(t<0.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    }
    
   
	uv=float3(u,v,0.0);
    n=normalize(dett);
}

void RDisk(float3 tp0, float3 dp1, float3 dp2, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
	float3 dp0=rp0-tp0;

	float3 dett =cross(dp1,dp2);
	float3 detuv=cross(dp0,rd);

	float det=(-1.0)/dot(dett,rd);

	float u=(dot(detuv,dp2))*det;
	float v=(dot(detuv,dp1))*det;
	t=(dot(dett ,dp0))*det;
    if(t<0.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    }
    
    if((u*u+v*v)>1.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    }
        
	uv=float3(u,v,0);    
    n=normalize(dett);
}

void RDDisk(float3 tp0, float3 np0, float r, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
	float3 dp0=rp0-tp0;

	float3 dp1;
	float3 dp2;
    np0=normalize(np0);

	if(abs(np0.x)<abs(np0.y))
		dp2=float3(1,0,0);
	else
		dp2=float3(0,1,0);
		
	dp1=normalize(cross(dp2,np0))*r;
	dp2=normalize(cross(dp1,np0))*r;
    
    
	float3 dett =cross(dp1,dp2);
	float3 detuv=cross(dp0,rd);

	float det=(-1.0)/dot(dett,rd);

	float u=(dot(detuv,dp2))*det;
	float v=(dot(detuv,dp1))*det;
	t=(dot(dett ,dp0))*det;
    if(t<0.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    }
    
    if((u*u+v*v)>1.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    }
        
	uv=float3(u,v,0);    
    n=normalize(dett);
}

void RCone(float3 p0, float r0, float3 p1, float r1, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)    
{
	float3 locX;
	float3 locY;
	float3 locZ=-(p1-p0)/(1.0-r1/r0);

    rp0-=p0-locZ;

	if(abs(locZ.x)<abs(locZ.y))
		locX=float3(1,0,0);
	else
		locX=float3(0,1,0);
		
	float len=length(locZ);
	locZ=normalize(locZ)/len;
	locY=normalize(cross(locX,locZ))/r0;
	locX=normalize(cross(locY,locZ))/r0;

	float3x3 tm;
	tm[0]=locX;
	tm[1]=locY;
	tm[2]=locZ;

    rd=rd*tm;	
    rp0=rp0*tm;
    	
	float dx=rd.x;
	float dy=rd.y;
	float dz=rd.z;

	float x0=rp0.x;
	float y0=rp0.y;
	float z0=rp0.z;

	float x02=x0*x0;
	float y02=y0*y0;
	float z02=z0*z0;

	float dx2=dx*dx;
	float dy2=dy*dy;
	float dz2=dz*dz;

	float det=(
		-2.0*x0*dx*z0*dz
        +2.0*x0*dx*y0*dy
        -2.0*z0*dz*y0*dy
        +dz2*x02
        +dz2*y02
        +dx2*z02
        +dy2*z02
        -dy2*x02
        -dx2*y02
        );
    

    if(det<0.0)
    {
		t=MAX_RAY_LENGTH;
        return;
    }

	float t0=(-x0*dx+z0*dz-y0*dy-sqrt(abs(det)))/(dx2-dz2+dy2);
	float t1=(-x0*dx+z0*dz-y0*dy+sqrt(abs(det)))/(dx2-dz2+dy2);

	t=t0;
	if(t<0.0)
    {
		t=MAX_RAY_LENGTH;
        return;
    }

	float3 pt=rp0+t*rd;

	if(pt.z>1.0)
    {
		t=MAX_RAY_LENGTH;
        return;
    }
        
    if(pt.z<r1/r0)
    {
		t=MAX_RAY_LENGTH;
        return;
    }

	n=float3(pt);
    uv.z=0.0;
    uv.y=n.z;
	n.z=0.0;
	n=normalize(n);
    uv.x=atan(n.x,n.y)/2.0/PI;
	n.z=-pt.z/abs(pt.z);
	n=normalize(n);
    n=tm*n;
    n=normalize(n);
}

void RSph(float3 p0, float r, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
	float3 l=p0-rp0;
	float tc=dot(l,rd);
	if(tc<0.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    };

    float d2=r*r+tc*tc-dot(l,l);

	if(d2<0.0)
    {
        t=MAX_RAY_LENGTH;
        return;
    };

	float thc=sqrt(d2);
    t=tc-thc;
    float3 p=rp0+rd*t;
    n=normalize(p-p0);
    uv.x=atan(n.x,n.z)/2.0/PI;
    uv.y=asin(n.y)/PI;
    uv.z=0.0;
}

void RCyl(float3 p0, float3 p1, float r, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
	float r2=r*r;


	float3 dp=p1-p0;
	float3 dpt=dp/dot(dp,dp);

	float3 ao=rp0-p0;
	float3 aoxab=cross(ao,dpt);
	float3 vxab=cross(rd,dpt);
	float ab2=dot(dpt,dpt);
	float a=2.0*dot(vxab,vxab);
	float ra=1.0/a;
	float b=2.0*dot(vxab,aoxab);
	float c=dot(aoxab,aoxab)-r2*ab2;

	float det=b*b-2.0*a*c;

	if(det<0.0)
    {
		t=MAX_RAY_LENGTH;
        return;
     }


	det=sqrt(det);

    float t0=(-b+det)*ra;
	float t1=(-b-det)*ra;

	if(t0>t1)
	{
		float temp=t1;
		t1=t0;
		t0=temp;
	}
	float d=t0;
	if(d<0.0)
    {
		t=MAX_RAY_LENGTH;
        return;
    }

	float3 ip=rp0+rd*d;
	float3 lp=ip-p0;
	float ct=dot(lp,dpt);
	if((ct<0.0)||(ct>1.0))
	{
		d=t1;
		if(d<0.0)
        {
            t=MAX_RAY_LENGTH;
            return;
        }

		ip=rp0+rd*d;
		float3 lp=ip-p0;
        float ct=dot(lp,dpt);
		if((ct<0.0)||(ct>1.0))
        {
        	t=MAX_RAY_LENGTH;
            return;
        }
	}

	t=d;
    n=normalize(ip-(p0+dp*ct));
    uv.y=ct;
	uv.x=n.x;
    uv.z=0.0;
}

void RRCone(float3 p0, float r0, float3 p1, float r1, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
 float3 l  = p1-p0;
 float ld = length(l);
 l=l/ld;
 float d=r0-r1;
 float sa = d/ld;
 float h0=r0*sa;
 float h1=r1*sa;
 float cr0 = sqrt(r0*r0-h0*h0);
 float cr1 = sqrt(r1*r1-h1*h1);
 float3 coneP0=p0+l*h0;
 float3 coneP1=p1+l*h1;
    
    float t0=MAX_RAY_LENGTH;
    {
        float t1;
        float3 uv1;
        float3 n1;
	    RCone(coneP0,cr0,coneP1,cr1,rp0,rd,t1,uv1,n1);
        if(t1<t0)
        {
            t0=t1;
            uv=uv1;
            n=n1;
        }
	    RSph(p0,r0,rp0,rd,t1,uv1,n1);
        if(t1<t0)
        {
            t0=t1;
            uv=uv1;
            n=n1;
        }
	    RSph(p1,r1,rp0,rd,t1,uv1,n1);
        if(t1<t0)
        {
            t0=t1;
            uv=uv1;
            n=n1;
        }
    }
    t=t0;
    
}

float3x3 Transpose(in float3x3 m)
{
	float3 i0 = m[0];
	float3 i1 = m[1];
	float3 i2 = m[2];
	float3x3 o=float3x3(
                 float3(i0.x, i1.x, i2.x),
                 float3(i0.y, i1.y, i2.y),
                 float3(i0.z, i1.z, i2.z)
                 );
	return o;
}
void REll(float3 p0, float3 r0, float3 r1, float3 r2, float3 rp0, float3 rd, out float t, out float3 uv, out float3 n)
{
    float3 irp0=rp0-p0;
//	float3 ir0=r0;
//	float3 ir1=r1;
//	float3 ir2=r2;

    float3 ir0=r0/dot(r0,r0);
	float3 ir1=r1/dot(r1,r1);
	float3 ir2=r2/dot(r2,r2);
//	r0=normalize(r0)/length(r0);
//	r1=normalize(r1)/length(r1);
//	r2=normalize(r2)/length(r2);

	float3x3 tm;
	tm[0]=ir0;
	tm[1]=ir1;
	tm[2]=ir2;

//    tm=Transpose(tm);
    
    float3 ird=rd*tm;	
    irp0=irp0*tm;	
    float t1=MAX_RAY_LENGTH;
    float3 uv1;
    float3 n1;
    float lr=length(ird);
    ird=normalize(ird);
    RSph(float3(0.0,0.0,0.0),1.0,irp0,ird,t1,uv1,n1);
    n=normalize(tm*n1);
    t=t1/lr;
    uv=uv1;
}

void trace(float3 rp0, float3 rd, out float t, out float3 col, out float3 n)
{
    float t1=1000.0;
    float3 col1;
    float3 n1;

    {
    	RP(float3(0.0,-1.0,0.0),float3(-1.0,0.0,0.0),float3(0.0,0,1.0),rp0, rd, t1, col1, n1);
        float3 p=rp0+rd*t1;
    	col1=float3(floor(mod(floor(p.z), 2.0)));
        if(mod(floor(p.x),2.0)==0.0)
        {
            p/=2.0;
    		col1=float3(floor(mod(floor(p.x+p.z+0.25)+floor(p.z-p.x+0.25), 2.0)));
        }
            
    }

    t=t1;
    col=col1;
    n=n1;


    float3 coneP0=float3(0.0,0.0,0.0);
    float3 coneP1=float3(0.0,0.0,3.0);
    {
        float t1;
        float3 col1;
        float3 n1;
	    RCone(coneP0,2.0,coneP1,1.0,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            float x=mod(floor(col1.x*8.0)+floor(col1.y*4.0),2.0);
            col=float3(x,x,x);
            t=t1;
            n=n1;
        }
    }
    
    {
        float t2=MAX_RAY_LENGTH;
        float3 col2;
        float3 n2;
//	    RSph(float3(-1.0,1.0,sin(iGlobalTime)*4.0),1.0,rp0,rd,t1,col1,n1);
        REll(
            //float3(-1.0,1.0,sin(iGlobalTime)*4.0),
            float3(0.0,3.0,0.0),
             float3(1.0,0.0,0.0),
             float3(0.0,2.0,0.0),
             float3(0.0,0.0,1.0),
             rp0,rd,t2,col2,n2);
        if(t2<t)
        {
            float x=mod(floor(col2.x*8.0)+floor(col2.y*4.0),2.0);
            col=float3(x,x,x);
            t=t2;
            n=n2;
        }
    }
    
//return;
    
    {
        float t1;
        float3 col1;
        float3 n1;
	    RDDisk(coneP0,coneP1-coneP0,2.0,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }

    {
        float t1;
        float3 col1;
        float3 n1;
	    RDDisk(coneP1,coneP0-coneP1,1.0,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }
    
    

    {
        float t1;
        float3 col1;
        float3 n1;
	    RSph(float3(1.0,3.0,3.0),0.4,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }

    {
        float t1;
        float3 col1;
        float3 n1;
	    RDisk(float3(0.0,0.0,1.0),float3(0.0,-1.0,0.0),float3(1.0,0.0,0.0),rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }
    
    {
        float t1;
        float3 col1;
        float3 n1;
	    RCyl(float3(0.0,0.0,1.0),float3(0.0,0.0,-1.0),1.0,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }

    {
        float t1;
        float3 col1;
        float3 n1;
	    RRCone(float3(3.0,0.0,2.0),1.0,float3(3.0,0.0,1.0),0.5,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(0.0,1.0,0.0);
            t=t1;
            n=n1;
        }
    }
    
    {
        float t1;
        float3 col1;
        float3 n1;
	    RSph(float3(0.0,0.0,-1.0),1.0,rp0,rd,t1,col1,n1);
        if(t1<t)
        {
            col=float3(1.0,0.5,0.0);
            t=t1;
            n=n1;
        }
    }
}


void lit(in float3 p, in float3 rd, in float3 n, in float3 icol, out float3 col)
{
    float3 tolight=normalize(light-p);

    float diffuse=clamp(dot(tolight,n),0.0,1.0);
    
    float3 halfNormal=normalize(tolight-rd);

    float3 nr=n*dot(n,-rd);
    float3 refl=normalize(-rd+(nr+rd)*2.0);
    
    float fresnel=(1.0-dot(-rd,n));
    float RF=0.2;
    fresnel=RF+(1.0-RF)*pow(1.0-dot(-rd,n),5.0);
    diffuse*=1.0-fresnel;
    
    float spec1=clamp(dot(n,halfNormal),0.0,1.0);
    float spec2=clamp(dot(tolight,refl),0.0,1.0);
    
    spec1=pow(spec1,20.0);
    spec2=pow(spec2,80.0)*2.0;
    float spec=spec1+(1.0-spec1)*spec2;
    
    diffuse=pow(diffuse,1.5);

    float shadow=1.0;
    float t1;
    float3 cols;
    float3 ns;
    trace(p+tolight*0.001,tolight,t1,cols,ns);
    if(t1<100.0)
    {
       shadow=0.0;
       spec=0.0;
    }
    diffuse*=shadow;
    
    col=icol;
    
    col*=(0.2+diffuse*0.8);
    col=clamp(col+(0.5+col*0.5)*spec1*(0.2+fresnel),0.0,1.0);
    col=mix(col,float3(1.0,1.0,1.0), clamp(spec2*diffuse*(1.0+fresnel),0.0,1.0));
}

void shade(float3 rp0, float3 rd, out float t, out float3 col, out float3 n)
{
    trace(rp0,rd,t,col,n);
//	col=n*0.5+0.5;
    float3 tolight=normalize(light-(rp0+rd*t));
    float diffuse=clamp(dot(tolight,n),0.0,1.0);
    
    float3 halfNormal=normalize(tolight-rd);

    float3 nr=n*dot(n,-rd);
    float3 refl=normalize(-rd+(nr+rd)*2.0);
    
    float fresnel=(1.0-dot(-rd,n));
    float RF=0.2;
    fresnel=RF+(1.0-RF)*pow(1.0-dot(-rd,n),5.0);
    diffuse*=1.0-fresnel;
    
    float spec1=clamp(dot(n,halfNormal),0.0,1.0);
    float spec2=clamp(dot(tolight,refl),0.0,1.0);
    
    spec1=pow(spec1,20.0);
    spec2=pow(spec2,80.0)*2.0;
    float spec=spec1+(1.0-spec1)*spec2;
//	spec=spec*1.2;
    float3 pos;
//    pos=n;col=fract(pos*0.5+0.5);
//    pos=rp0+rd*t;col=fract(pos*4.0);
//    return;
    
	float shadow=1.0;
//    if(false)
    {
	    float t1;
    	float3 col1;
    	float3 n1;
        float3 pos=rp0+t*rd+n*0.001;
        trace(pos,normalize(light-pos), t1, col1, n1);
        if(t1<100.0)
        {
            shadow=0.0;
            spec=0.0;
        }
    }
    spec1*=shadow;
    spec2*=shadow;
    diffuse=pow(diffuse,1.5);
    diffuse*=shadow;
    col*=(0.2+diffuse*0.8);
//	return;
//    if(false)
    {
        
	    float t1;
    	float3 col1;
    	float3 n1;
        float3 pos=rp0+t*rd+n*0.001;
        trace(pos,refl, t1, col1, n1);
        float3 col2=col1;
        lit(pos+refl*t1, refl, n1, col2, col1);
	    
        float3 fogcol=mix(float3(0.87,0.8,0.83),float3(0.3,0.6,1.0),1.0-(1.0-refl.y)*(1.0-refl.y));
    	col1=mix(fogcol,col1,clamp(1.8/(exp(t1*0.25)),0.0,1.0));
        
//        if(t1<100.0)
        {
//            col+=(col*0.5+0.5)*col1/exp(t1*0.05)*clamp(dot(tolight,n1),0.0,1.0)*fresnel;
            col+=col1*(0.3+fresnel*0.7);
        }
    }

    col=clamp(col+(0.5+col*0.5)*spec1*(0.2+fresnel),0.0,1.0);
    col=mix(col,float3(1.0,1.0,1.0), clamp(spec2*diffuse*(1.0+fresnel),0.0,1.0));

    float3 fogcol=mix(float3(0.87,0.8,0.83),float3(0.3,0.6,1.0),1.0-(1.0-rd.y)*(1.0-rd.y));
    col=mix(fogcol,col,clamp(1.8/(exp(t*0.025)),0.0,1.0));
    
    
//    col*=(0.5+shadow*0.5);
}

void main(void)
{
    float T=iGlobalTime*0.45;
    
    light.x=cos(T)*10.0;
    light.z=sin(T)*10.0;
    light.y=5.0;
    
    float mposx=iMouse.x;
    float mposy=iMouse.y;
    if(iMouse.z<0.0)mposx=-iMouse.z;
    if(iMouse.w<0.0)mposy=-iMouse.w;
    
    float a1=-(mposy/iResolution.y)*PI/2.1+0.1;
    float a2=mposx/iResolution.x*PI*2.0-0.3;
    campos.y=sin(a1)*campos.x;
    float camx=cos(a1)*6.0;
    campos.x=cos(a2)*camx;
    campos.z=sin(a2)*camx;
    campos+=look_at;
    
    forward=normalize(look_at-campos);
    right=normalize(cross(up,forward));
    up=normalize(cross(forward,right));
    
	float2 scr = gl_FragCoord.xy /iResolution.xy;
    scr=2.0*scr-1.0;
    float2 scruv=scr;
    scr.x*=(iResolution.x/iResolution.y);
    
    float3 ray=pixelRay(scr);

    if (!setup_ray(eye, dir, campos, ray)) return;

    float t;
    float3 col=float3(0.0,0.0,0.0);
    float3 n;

    
//    shade(campos, ray, t, col, n);
    float w=0.0;
    float3 col1;
    const float nx=2.0;
    const float ny=2.0;
                      
	for(float i=0.0;i<ny;i+=1.001)
    {
		for(float j=0.0;j<nx;j+=1.001)
        {
	    //	shade(campos, pixelRay(scr+float2((j*2.0/(nx+0.0))/iResolution.y,(i*2.0/(ny+0.0))/iResolution.y)), t, col1, n);
	    	shade(campos, ray, t, col1, n);
            w=w+1.0;
            col+=col1;
        }
    }
  	col=col/w;
    col=col-0.25*dot(scruv.xy*abs(scruv.xy),scruv.xy);
    //gl_FragColor = float4(col,1.0);
    write_pixel(dir, 0., col);
}
