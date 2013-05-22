// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//
// I can't recall where I learnt about this fractal.
//
// Coloring and fake occlusions are done by orbit trapping, as usual (unless somebody has invented
// something new in the last 4 years that i'm unaware of, that is)
//
// From: https://www.shadertoy.com/view/4ds3zn

// Camera position, direction and eps from vertex shader.
varying vec3 eye;
varying vec3 dir;
varying float zoom;

uniform float xres, yres, time, speed;
uniform int iters;    // {min=1 max=1000} Number of fractal iterations.

vec2      iResolution = vec2(xres, yres);     // viewport resolution (in pixels)
float     iGlobalTime = time;     // shader playback time (in seconds)

uniform vec3 par[1];
#define ss par[0].x  // {min=.1 max=2.0 step=.001}

uniform float focus;  // {min=-10 max=30 step=.01} Focal plane devation from 30x speed.
bool setup_stereo(inout vec3 eye_in, inout vec3 dp) {
#if !defined(ST_NONE)
#if defined ST_OCULUS
  float halfx = xres / 2.0;

  vec2 q;
  if (sign(speed) < 0.0) {
    // left. 45 pixel shift towards center. Eyeballed.
    q = (gl_FragCoord.xy - vec2(focus + 45.0, 0.0)) / vec2(halfx, yres);
  } else {
    // right. 45 pixel shift towards center.
    q = (gl_FragCoord.xy - vec2(halfx - focus - 45.0, 0.0)) / vec2(halfx, yres);
  }
  vec2 p = -1.0 + 2.0 * q;

  // Oculus barrel distort parameters.
  vec3 oculus_warp = vec3(1.0, 0.22, 0.24);  // k0, k1, k2
  vec2 oculus_scale = vec2(0.3, 0.35);  // x/y ratio eyeballed
  float r2 = dot(p, p);  // Radius squared, from center.
  p *= oculus_scale * dot(oculus_warp, vec3(1.0, r2, r2*r2));
  if (dot(p, p) > 0.10) { 
    //discard;  // Don't waste time on pixels we can't see.
    return false;
  }

  // Shift eye position, abs(speed) is half inter-occular distance.
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0.0, 0.0, 0.0));
  eye_in = eye + eye_d;

  // Note: no asymmetric frustum for Rift.
  dp = normalize(vec3(gl_ModelViewMatrix * vec4(p, 0.35, 0.0)));  // z value determines fov. Eyeballed.
#else
#if defined(ST_INTERLACED)
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * abs(speed), 0, 0, 0));
#else
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4(speed, 0, 0, 0));
#endif
  eye_in = eye + eye_d;
  // Construct asymmetric frustum.
  dp = normalize(dir * (focus + 30.0) * abs(speed) - eye_d);
#endif // ST_OCULUS
#else  // ST_NONE
  eye_in = eye;
  dp = normalize(dir);
#endif
  return true;
}

vec4 orb; 
//float ss = 1.1;
float map( vec3 p )
{
	float scale = 1.0;

	orb = vec4(1000.0); 
	
	for( int i=0; i<iters;i++ )
	{
		p = -1.0 + 2.0*fract(0.5*p+0.5);

		float r2 = dot(p,p);
		
		orb = min( orb, vec4(abs(p),r2) );
		
		float k = max(ss/r2,0.1);
		p     *= k;
		scale *= k;
	}
	
	return 0.25*abs(p.y)/scale;
}

float trace( in vec3 ro, in vec3 rd )
{
	float maxd = 100.0;
	float m_zoom = zoom * .5 / xres;
	float precis = 0.001;
  float h=precis*2.0;
	float t = 0.0;
  for( int i=0; i<200; i++ )
    {
        if( abs(h)<precis||t>maxd ) break;
        t += h;
	      precis = m_zoom * t;
	      h = map( ro+rd*t );
    }

  if( t>maxd ) t=-1.0;
  return t;
}

vec3 calcNormal( in vec3 pos )
{
	vec3  eps = vec3(.00001,0.0,0.0);
	vec3 nor;
	nor.x = map(pos+eps.xyy) - map(pos-eps.xyy);
	nor.y = map(pos+eps.yxy) - map(pos-eps.yxy);
	nor.z = map(pos+eps.yyx) - map(pos-eps.yyx);
	return normalize(nor);
}

void main(void)
{
#if 0
  vec2 p = -1.0 + 2.0*gl_FragCoord.xy / iResolution.xy;
  p.x *= iResolution.x/iResolution.y;

	float time = iGlobalTime*0.25 + 0.01*iMouse.x;
  ss = 1.1 + 0.5*smoothstep( -0.3, 0.3, cos(0.1*iGlobalTime) );
	
  // camera
	vec3 ro = vec3( 2.8*cos(0.1+.33*time), 0.4 + 0.30*cos(0.37*time), 2.8*cos(0.5+0.35*time) );
	vec3 ta = vec3( 1.9*cos(1.2+.41*time), 0.4 + 0.10*cos(0.27*time), 1.9*cos(2.0+0.38*time) );
	float roll = 0.2*cos(0.1*time);
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(roll), cos(roll),0.0);
	vec3 cu = normalize(cross(cw,cp));
	vec3 cv = normalize(cross(cu,cw));
	vec3 rd = normalize( p.x*cu + p.y*cv + 2.0*cw );
#else
	vec3 ro;
	vec3 rd;
	if (!setup_stereo(ro, rd)) {
		gl_FragColor = vec4(0.0);
		return;
	}
#endif

  // trace	
	vec3 col = vec3(0.0);
	float t = trace( ro, rd );
	if( t>0.0 )
	{
		vec4 tra = orb;
		vec3 pos = ro + t*rd;
		vec3 nor = calcNormal( pos );
		
		// lighting
					vec3  light1 = vec3(  0.577, 0.577, -0.577 );
					vec3  light2 = vec3( -0.707, 0.000,  0.707 );
		float key = clamp( dot( light1, nor ), 0.0, 1.0 );
		float bac = clamp( 0.2 + 0.8*dot( light2, nor ), 0.0, 1.0 );
		float amb = (0.7+0.3*nor.y);
		float ao = pow( clamp(tra.w*2.0,0.0,1.0), 1.2 );

		vec3 brdf  = 1.0*vec3(0.40,0.40,0.40)*amb*ao;
		brdf += 1.0*vec3(1.00,1.00,1.00)*key*ao;
		brdf += 1.0*vec3(0.40,0.40,0.40)*bac*ao;

					// material		
		vec3 rgb = vec3(1.0);
		rgb = mix( rgb, vec3(1.0,0.80,0.2), clamp(6.0*tra.y,0.0,1.0) );
		rgb = mix( rgb, vec3(1.0,0.55,0.0), pow(clamp(1.0-2.0*tra.z,0.0,1.0),8.0) );

		// color
		col = rgb*brdf*exp(-0.2*t);
	}

	col = sqrt(col);
	
	col = mix( col, smoothstep( 0.0, 1.0, col ), 0.25 );
	
	gl_FragColor=vec4(col,1.0);
}
