// from: https://www.shadertoy.com/view/4dB3DW

#include "setup.inc"
#line 5

float iGlobalTime = 0.0;
vec2 iResolution = vec2(xres, yres);

#define maxSteps 76.0
#define treshold 0.001
#define maxdist 20.0
#define shadowsteps 10.0
#define aosteps 10.0
#define pi acos(-1.)

float perlin(vec3 p) {
	vec3 i = floor(p);
	vec4 a = dot(i, vec3(1., 57., 21.)) + vec4(0., 57., 21., 78.);
	vec3 f = cos((p-i)*pi)*(-.5)+.5;
	a = mix(sin(cos(a)*a),sin(cos(1.+a)*(1.+a)), f.x);
	a.xy = mix(a.xz, a.yw, f.y);
	return mix(a.x, a.y, f.z);
	}


vec2 map(vec3 p ) {
	p=p*0.5 + vec3(-1.5,-1.5+cos(iGlobalTime)*0.5,-iGlobalTime);
	return vec2( mix(sin(p.x),sin(p.y),sin(p.z)) - perlin(p) , 1.0);
	}



vec2 rot(vec2 k, float t) {
	return vec2(cos(t)*k.x-sin(t)*k.y,sin(t)*k.x+cos(t)*k.y);
	}

vec3 cNor(vec3 p ) {
	vec3 e=vec3(0.001,0.0,0.0);
	return normalize(vec3( map(p+e.xyy).x - map(p-e.xyy).x, map(p+e.yxy).x - map(p-e.yxy).x, map(p+e.yyx).x - map(p-e.yyx).x ));
	}

float cShd(vec3 ro, vec3 rd, float k ) {
	float res = 1.0;
	for(float i=1.0; i<shadowsteps; i+=1.0){
		float f=shadowsteps/i;
        float h = map(ro + rd*f).x;
        if( h<0.001 ) { res=0.0; break; }
        res = min( res, k*h/f );
    }
    return res;
}

float calcAO(vec3 pos, vec3 nor ){
	float totao = 0.0;
    float sca = 1.0;
    for( float aoi=1.0; aoi<aosteps; aoi+=1.0 ) {
        float hr = 0.05*aoi;
        float dd = map( pos + nor*hr ).x;
        totao += -(dd-hr)*sca;
        sca *= 0.75;
    }
    return clamp( 1.0 - 4.0*totao, 0.0, 1.0 );
}


vec3 tmap(vec3 p, vec3 n) {
	float f=abs(perlin(n*6.0))*0.2+0.2;
	vec3 col=f*vec3(0.8)+0.5;
	return col*vec3(0.6,0.3,0.2)+normalize(p)*0.2;	
	}


void main(void)	{
	vec2 ps=(gl_FragCoord.xy/iResolution.xy);
	vec3 rd=normalize( vec3( (-1.0+2.0*ps)*vec2(1.0,1.0), 1.0));
	vec3 ro=vec3(0.0, 0.0, -0.5);
	vec3 lig=vec3(0.0,0.0,-0.5);

#if 0
	vec4 mouse=iMouse*0.01;
	lig.xz=rot(lig.xz, mouse.x);
	lig.xy=rot(lig.xy, mouse.y);
	ro.xz=rot(ro.xz, mouse.x);
	ro.xy=rot(ro.xy, mouse.y);
	rd.xz=rot(rd.xz, mouse.x);
	rd.xy=rot(rd.xy, mouse.y);	
#else
  if (!setup_ray(eye, dir, ro, rd)) return;
#endif
	
	//march
	float f=0.0;
	vec2 t=vec2(treshold,f);
	for(float i=0.0; i<1.0; i+=1.0/maxSteps){
        t= map(ro + rd*t.x);
		f+=t.x;
		t.x=f;
        if( abs(t.x)<treshold || t.x>maxdist ) break; 
		}
	if (t.x>maxdist) t.y=0.0;

	//draw
	vec3 col = vec3(0.0);
	if (t.y>0.5) {
		
		lig=normalize(lig);
		vec3 pos = ro + rd*t.x;
		vec3 nor = cNor(pos);
		float ao = calcAO( pos, nor );

		
		float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0 );					
		float dif = clamp( dot( nor, lig ), 0.0, 1.0 );				
		float bac = clamp( dot( nor, vec3(-lig.x,lig.y,-lig.z)), 0.0, 1.0 );		//1.0

		float sh = cShd( pos, lig, 7.0 );	

		col = 0.20*amb*vec3(0.10,0.10,0.10)*ao;						
		col += 0.20*bac*vec3(0.15,0.15,0.15)*ao;					
		col += 1.00*dif*vec3(0.80,0.80,0.80);					

		float spe = sh*pow(clamp( dot( lig, reflect(rd,nor) ), 0.0, 1.0 ) ,16.0 );	//1.0
		float rim = ao*pow(clamp( 1.0+dot(nor,rd),0.0,1.0), 2.0 );	
	vec3 oc=tmap(pos,nor);
		
		col =oc*col + vec3(1.0)*col*spe + 0.3*rim*(0.5+0.5*col);	
		col*=exp(.07*f); col*=2.;
	
	} 
		
#if 0
	gl_FragColor=vec4( col, 1.0);
#else
  write_pixel(dir, t.x, col);
#endif
	}
