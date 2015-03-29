// Space Jewels. December 2014
// https://www.shadertoy.com/view/llX3zr

#include "setup.inc"
#line 6

uniform int color_iters, iters;

vec2 iResolution = vec2(xres, yres);
#define iGlobalTime time
vec4 iMouse = vec4(0.);

//--------------------------------------------------------------------------
#define PI 3.141596
vec3 sunLight  = normalize( vec3(  0.5, 0.3,  0.3 ) );
const vec3 sunColour = vec3(1., .9, .85);
#define CSize  vec3(.808, .8, 1.137)
const vec3 fogColour = vec3(0.09, 0.08, 0.07);
float gTime;

//----------------------------------------------------------------------------------------
float Hash(vec2 p)
{
	p  = fract(p * vec2(5.3983, 5.4427));
    p += dot(p.yx, p.xy + vec2(21.5351, 14.3137));
	return fract(p.x * p.y * 95.4337);
}

//----------------------------------------------------------------------------------------
vec3 Colour( vec3 p)
{
	float col	= 0.0;
    float r2	= dot(p,p);
		
	for( int i=0; i < color_iters;i++ )
	{
		vec3 p1= 2.0 * clamp(p, -CSize, CSize)-p;
		col += abs(p.z-p1.z);
		p = p1;
		r2 = dot(p,p);
		float k = max((1.1)/(r2), .03);
		p *= k;
	}
	return (0.5+0.5*sin(col*vec3(1.647,-1.0,4.9)));
}

//--------------------------------------------------------------------------

float Map( vec3 p )
{
	float scale = 1.0;
	
	for( int i=0; i < iters;i++ )
	{
		p = 2.0*clamp(p, -CSize, CSize) - p;
		float r2 = dot(p,p);
		float k = max((1.1)/(r2), .03);
		p     *= k;
		scale *= k;
	}
	float l = length(p.xy);
	float rxy = l - float(iters) - 1.0;
	float n = l * p.z;
	rxy = max(rxy, -(n) / (length(p))-.1);
	return (rxy) / abs(scale);
}

//--------------------------------------------------------------------------
float SphereRadius(float t)
{
	if (t< 0.4) t=  abs(t-0.4) * 5.5;
	t = t*0.05;
	return clamp(t*t, 5.0/iResolution.x, 400.0/iResolution.x);
}

//--------------------------------------------------------------------------
float Shadow( in vec3 ro, in vec3 rd)
{
	float res = 1.0;
    float t = 0.2;
	float h;
	
    for (int i = 0; i < 5; i++)
	{
		h = Map( ro + rd*t );
		res = min(7.0*h / t, res);
		t += h*.5+.01;
	}
    return max(res, 0.0);
}

//--------------------------------------------------------------------------
vec3 DoLighting(in vec3 mat, in vec3 pos, in vec3 normal, in vec3 eyeDir, in float d)
{
	float sh = Shadow(pos+normal*.01, sunLight);
    // Light surface with 'sun'...
	vec3 col = mat * sunColour*(max(dot(sunLight,normal), 0.0)) *sh;
    
    normal = reflect(eyeDir, normal); // Specular...
    col += pow(max(dot(sunLight, normal), 0.0), 25.0)  * sunColour * 1.5 *sh;
    // Ambient..
    col += mat * .1 * max(normal.z, 0.0);
    col = mix(fogColour,col, min(exp(-d*d*.05), 1.0));
    
	return col;
}


//--------------------------------------------------------------------------
vec3 GetNormal(vec3 p, float sphereR)
{
	vec2 eps = vec2(sphereR, 0.0);
	return normalize( vec3(
           Map(p+eps.xyy) - Map(p-eps.xyy),
           Map(p+eps.yxy) - Map(p-eps.yxy),
           Map(p+eps.yyx) - Map(p-eps.yyx) ) );
}

//--------------------------------------------------------------------------
float Scene(in vec3 rO, in vec3 rD, inout vec4 aStack1, inout vec4 aStack2, inout vec4 dStack1,  inout vec4 dStack2)
{
    //float t = 0.0;
	float t = .1 * Hash(gl_FragCoord.xy*fract(iGlobalTime));
	float  alphaAcc = 0.0;
	vec3 p = vec3(0.0);

	for( int j=0; j < 80; j++ )
	{
		if (alphaAcc >= 1. || t > 15.0) break;
		p = rO + t*rD;
		float sphereR = SphereRadius(t);
		float h =Map(p);
		if( h < sphereR)
		{
			// Accumulate the alphas...
			float alpha = clamp((1.0 - alphaAcc) * ((sphereR-h) / sphereR), 0.0, 1.0);
            // If high enough to contribute nicely...
            if (alpha >= (1.0 / 8.0))
            {
                // put it on the 2 stacks, alpha and depth...
                aStack2.yzw = aStack2.xyz; aStack2.x = aStack1.w;
                aStack1.yzw = aStack1.xyz; aStack1.x = alpha;
                dStack2.yzw = dStack2.xyz; dStack2.x = dStack1.w;
                dStack1.yzw = dStack1.xyz; dStack1.x = t;
				alphaAcc += alpha;	
            }
            
		}
		t +=  h*.8 + t*0.001;
        
	}
    if (t > 15.0)
    {
        // Make sure the far distance is hit properly...
		aStack2.yzw = aStack2.xyz; aStack2.x = aStack1.w;
		aStack1.yzw = aStack1.xyz; aStack1.x = (1.0 - alphaAcc);
		dStack2.yzw = dStack2.xyz; dStack2.x = dStack1.w;
		dStack1.yzw = dStack1.xyz; dStack1.x = t;
		alphaAcc += aStack1.x;
    }
    alphaAcc = clamp(alphaAcc, 0.0, 1.0);
    
	return alphaAcc;
}


//--------------------------------------------------------------------------
vec3 PostEffects(vec3 rgb, vec2 xy)
{
	// Gamma first...
	rgb = pow(rgb, vec3(0.45));

	// Then...
	#define CONTRAST 1.3
	#define SATURATION 1.3
	#define BRIGHTNESS 1.2
	rgb = mix(vec3(.5), mix(vec3(dot(vec3(.2125, .7154, .0721), rgb*BRIGHTNESS)), rgb*BRIGHTNESS, SATURATION), CONTRAST);

	// Vignette...
	rgb *= .5+0.5*pow(180.0*xy.x*xy.y*(1.0-xy.x)*(1.0-xy.y), 0.3 );	

	return clamp(rgb, 0.0, 1.0);
}

//--------------------------------------------------------------------------
vec3 TexCube( sampler2D sam, in vec3 p, in vec3 n )
{
	vec3 x = texture2D( sam, p.yz ).xyz;
	vec3 y = texture2D( sam, p.zx ).xyz;
	vec3 z = texture2D( sam, p.xy ).xyz;
	return (x*abs(n.x) + y*abs(n.y) + z*abs(n.z))/(abs(n.x)+abs(n.y)+abs(n.z));
}

//--------------------------------------------------------------------------
vec3 Albedo(vec3 pos, vec3 nor)
{
    vec3 col = TexCube(iChannel0, pos*.1, nor).xyz;
    col *= Colour(pos);
    return col;
}


//--------------------------------------------------------------------------
vec3 CameraPath( float t )
{
    vec3 p = vec3(-13.0 +3.4 * sin(t),-0.+4.5 * cos(t),-1.1+.4 * cos(2.3*t+2.0) );
	return p;
} 
    

//--------------------------------------------------------------------------
void main(void)
{
	float m = (iMouse.x/iResolution.x)*20.0;
	float gTime = ((iGlobalTime+23.)*.2+m);
  vec2 xy = gl_FragCoord.xy / iResolution.xy;
	vec2 uv = (-1.0 + 2.0 * xy) * vec2(iResolution.x/iResolution.y,1.0);
	
	float hTime = mod(gTime+1.95, 2.0);
	
	vec3 cameraPos 	= CameraPath(gTime + 0.0);
	vec3 camTarget 	= vec3 (-12., -.0, -2.0);

  float roll = 0.0;//clamp(GetAngle(v1 , v2), -.8, .8);
    

	vec3 cw = normalize(camTarget-cameraPos);
	vec3 cp = vec3(sin(roll), cos(roll),180.0);
	vec3 cu = normalize(cross(cw,cp));
	vec3 cv = cross(cu,cw);
	vec3 ddir = normalize(uv.x*cu + uv.y*cv + 1.1*cw);

	vec3 col = vec3(0.0);
	
  vec4 dStack1 = vec4(-200.0);
	vec4 dStack2 = vec4(-200.0);
  vec4 aStack1 = vec4(.0);
  vec4 aStack2 = vec4(.0);

  if (!setup_ray(eye, dir, cameraPos, ddir)) return;

	float alpha = Scene(cameraPos, ddir, aStack1, aStack2, dStack1, dStack2);
    
    // Render both stacks...
    for (int i = 0; i < 4; i++)
    {
        float d = dStack1[i];
        
        if (d < .0) continue;
        float sphereR = SphereRadius(d);
        vec3 pos = cameraPos + ddir * d;
    		vec3 normal = GetNormal(pos, sphereR);
        vec3 c = Albedo(pos, normal);
        col += DoLighting(c, pos, normal, ddir, d) * aStack1[i];
    }
 
    for (int i = 0; i < 4; i++)
    {
        float d = dStack2[i];
	    	if (d < .0) continue;
    	
        float sphereR = SphereRadius(d);
        vec3 pos = cameraPos + ddir * d;
		    vec3 normal = GetNormal(pos, sphereR);
		    vec3 c = Albedo(pos, normal);
        col += DoLighting(c, pos, normal, ddir, d) * aStack2[i];
    }
   
	  col = PostEffects(col, xy) * smoothstep(.0, 2.0, iGlobalTime);	
	
    write_pixel(dir, 0., col);
	//gl_FragColor=vec4(col,1.0);
}

//--------------------------------------------------------------------------
