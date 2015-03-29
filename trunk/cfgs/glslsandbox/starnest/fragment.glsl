// from: http://glsl.heroku.com/e#19539.0

#include "setup.inc"
#line 4

uniform int iters;    // {min=1 max=1000} Number of fractal iterations.
uniform int color_iters;    // {min=1 max=1000} Number of fractal iterations.

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
vec2 mouse = vec2(0., 0.);
vec2 resolution = vec2(xres, yres);

// Star Nest by Pablo Román Andrioli
// Modified a lot.

// This content is under the MIT License.

#define iterations color_iters
#define formuparam 0.530

#define volsteps iters
#define stepsize 0.2

#define zoom   0.800
#define tile   0.850
#define speed  0.01

#define brightness 0.0015
#define darkmatter 0.400
#define distfading 0.760
#define saturation 0.800


void main(void)
{
	//get coords and direction
	vec2 uv=gl_FragCoord.xy/resolution.xy-.5;
	uv.y*=resolution.y/resolution.x;
	vec3 ddir=vec3(uv*zoom,1.);
	
	float a2=speed+.5;
	float a1=0.0;
	mat2 rot1=mat2(cos(a1),sin(a1),-sin(a1),cos(a1));
	mat2 rot2=rot1;//mat2(cos(a2),sin(a2),-sin(a2),cos(a2));
	ddir.xz*=rot1;
	ddir.xy*=rot2;
	
	//from.x-=time;
	//mouse movement
	vec3 from=vec3(0.,0.,0.);
	from+=vec3(.05,.05,-2.);
	
	//from.z+=time*.01;
	
	from.x-=mouse.x;
	from.y-=mouse.y;
	
	from.xz*=rot1;
	from.xy*=rot2;

  if (!setup_ray(eye, dir, from, ddir)) return;
	
	//volumetric rendering
	float s=.4,fade=.2;
	vec3 v=vec3(0.4);
	for (int r=0; r<volsteps; r++) {
		vec3 p=from+s*ddir*.5;
		p = abs(vec3(tile)-mod(p,vec3(tile*2.))); // tiling fold
		float pa,a=pa=0.;
		for (int i=0; i<iterations; i++) { 
			p=abs(p)/dot(p,p)-formuparam; // the magic formula
			a+=abs(length(p)-pa); // absolute sum of average change
			pa=length(p);
		}
		float dm=max(0.,darkmatter-a*a*.001); //dark matter
		a*=a*a*2.; // add contrast
		if (r>3) fade*=1.-dm; // dark matter, don't render near
		//v+=vec3(dm,dm*.5,0.);
		v+=fade;
		v+=vec3(s,s*s,s*s*s*s)*a*brightness*fade; // coloring based on distance
		fade*=distfading; // distance fading
		s+=stepsize;
	}
	v=mix(vec3(length(v)),v,saturation); //color adjust

	//gl_FragColor = vec4(v*.01,1.);	
  write_pixel(dir, 1.0, v*.01);
}
