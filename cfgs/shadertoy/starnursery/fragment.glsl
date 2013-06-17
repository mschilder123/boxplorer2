// Built from the basics of'Clouds' Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

// Edited by Dave Hoskins into "Star Nursery"
// V.1.1 Some speed up in the ray-marching loop.

// boxplorer i/o
varying vec3 eye, dir;
uniform float xres, yres, speed, time;

// Map to shadertoy expected vars
float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);

// boxplorer 3d hackery
uniform float focus;  // {min=-10 max=30 step=.01} Focal plane devation from 30x speed.
bool setup_stereo(out vec3 eye_in, out vec3 dp) {
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
  vec3 eye_d = vec3(gl_ModelViewMatrix * vec4( 2.0 * (fract(gl_FragCoord.y * 0.5) - .5) * speed, 0, 0, 0));
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


mat3 m = mat3( 0.00,  0.90,  0.60,
              -0.90,  0.36, -0.48,
              -0.60, -0.48,  0.34 );
//float time = iGlobalTime+5.4;

//----------------------------------------------------------------------
float hash( float n )
{
    return fract(sin(n)*43758.5453123);
}

//----------------------------------------------------------------------
float noise( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);

    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*57.0;

    float res = mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                    mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y);

    return res;
}

//----------------------------------------------------------------------
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
	
    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*57.0 + 113.0*p.z;

    float res = mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
    return res;
}

//----------------------------------------------------------------------
float fbm( vec3 p )
{
    float f;
    f  = 1.600*noise( p ); p = m*p*2.02;
    f += 0.3500*noise( p ); p = m*p*2.33;
    f += 0.2250*noise( p ); p = m*p*2.01;
    f += 0.0825*noise( p ); p = m*p*2.01;
    return f;
}

//----------------------------------------------------------------------
vec4 map( in vec3 p )
{
	float d = 0.2 - p.y;

	float f= fbm( p*1.0 - vec3(.4,0.3,-0.3)*time);
	d += 4.0 * f;

	d = clamp( d, 0.0, 1.0 );
	
	vec4 res = vec4( d );
	res.w = pow(res.y, .1);

	res.xyz = mix( .7*vec3(1.0,0.4,0.2), vec3(0.2,0.0,0.2), res.y * 1.);
	res.xyz = res.xyz + pow(abs(.95-f), 26.0) * 1.85;
	return res;
}


//----------------------------------------------------------------------
vec3 sundir = vec3(1.0,0.4,0.0);
vec4 raymarch( in vec3 ro, in vec3 rd )
{
	vec4 sum = vec4(0, 0, 0, 0);

	float t = 0.0;
	vec3 pos = vec3(0.0, 0.0, 0.0);
	for(int i=0; i<60; i++)
	{
		if (sum.a > 0.8 || pos.y > 9.0 || pos.y < -2.0) continue;
		pos = ro + t*rd;	

		vec4 col = map( pos );
		
		// Accumulate the alpha with the colour...
		col.a *= 0.08;
		col.rgb *= col.a;

		sum = sum + col*(1.0 - sum.a);	


    	t += max(0.1,0.08*t);
	}
	sum.xyz /= (0.003+sum.w);

	return clamp( sum, 0.0, 1.0 );
}

//----------------------------------------------------------------------
void main(void)
{
	vec2 q = gl_FragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0*q;
    p.x *= iResolution.x/ iResolution.y;
    vec2 mo = vec2(0.0, 0.0); //(-1.0 + 2.0 + iMouse.xy) / iResolution.xy;
    
    // Camera code...
    vec3 ro = 5.6*normalize(vec3(cos(2.75-3.0*mo.x), .4-1.3*(mo.y-2.4), sin(2.75-2.0*mo.x)));
	vec3 ta = vec3(.0, 5.6, 2.4);
    vec3 ww = normalize( ta - ro);
    vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww ));
    vec3 vv = normalize(cross(ww,uu));
    vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

    // Use boxplorer camera
    if (!setup_stereo( ro, rd )) {
      gl_FragColor = vec4(0);
      return;
    }

	// Ray march into the clouds adding up colour...
    vec4 res = raymarch( ro, rd );
	

	float sun = clamp( dot(sundir,rd), 0.0, 2.0 );
	vec3 col = vec3(0.2,0.2,0.3);
	col += .4*vec3(.4,.2,0.67)*sun;
	col = clamp(col, 0.0, 1.0);
	col += 0.43*vec3(.4,0.4,0.2)*pow( sun, 21.0 );
	
	// Do the stars...
	float v = 1.0/( 2. * ( 1. + rd.z ) );
	vec2 xy = vec2(rd.y * v, rd.x * v);
    float s = noise(rd.xz*134.);
	s += noise(rd.xz*370.);
	s += noise(rd.xz*870.);
	s = pow(s,19.0) * 0.00000001 * max(rd.y, 0.0);
	if (s > 0.0)
	{
		vec3 backStars = vec3((1.0-sin(xy.x*20.0+time*13.0*rd.x+xy.y*30.0))*.5*s,s, s); 
		col += backStars;
	}

	// Mix in the clouds...
	col = mix( col, res.xyz, res.w*1.3);
	
	    
    gl_FragColor = vec4( col, 1.0 );
}
