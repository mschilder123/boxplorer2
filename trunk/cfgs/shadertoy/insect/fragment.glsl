// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// From https://www.shadertoy.com/view/Mss3zM

// Boxplorer inputs
varying vec3 eye, dir;
varying float zoom;
uniform float xres;
uniform float yres;
uniform float time;
uniform float speed;
uniform sampler2D bg_texture;

// Boxplorer only provices 1 input texture; use it.
#define iChannel0 bg_texture
#define iChannel1 bg_texture
#define iChannel2 bg_texture
float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);

float hash( float n )
{
  return fract(sin(n)*158.5453);
}

float noise( in float x )
{
  float p = floor(x);
  float f = fract(x);

  f = f*f*(3.0-2.0*f);

  return mix( hash(p+0.0), hash(p+1.0),f);
}


float noise( in vec2 x )
{
  vec2 p = floor(x);
  vec2 f = fract(x);

  f = f*f*(3.0-2.0*f);
  float a = texture2D(iChannel1,(p+vec2(0.5,0.5))/64.0,-32.0).x;
  float b = texture2D(iChannel1,(p+vec2(1.5,0.5))/64.0,-32.0).x;
  float c = texture2D(iChannel1,(p+vec2(0.5,1.5))/64.0,-32.0).x;
  float d = texture2D(iChannel1,(p+vec2(1.5,1.5))/64.0,-32.0).x;
  float res = mix(mix( a, b,f.x), mix( c, d,f.x),f.y);

  return 2.0*res;
}

vec3 texturize( sampler2D sa, vec3 p, vec3 n )
{
  vec3 x = texture2D( sa, p.yz ).xyz;
  vec3 y = texture2D( sa, p.zx ).xyz;
  vec3 z = texture2D( sa, p.xy ).xyz;

  return x*abs(n.x) + y*abs(n.y) + z*abs(n.z);
}

//----------------------------------------------------------------

float terrainSoft( vec2 x )
{
  x += 100.0;
  x *= 0.6;

  vec2 p = floor(x);
  vec2 f = fract(x);

  f = f*f*(3.0-2.0*f);
  float a = texture2D(iChannel0,0.0+(p+vec2(0.5,0.5))/1024.0,-32.0).x;
  float b = texture2D(iChannel0,0.0+(p+vec2(1.5,0.5))/1024.0,-32.0).x;
  float c = texture2D(iChannel0,0.0+(p+vec2(0.5,1.5))/1024.0,-32.0).x;
  float d = texture2D(iChannel0,0.0+(p+vec2(1.5,1.5))/1024.0,-32.0).x;
  float r = mix(mix( a, b,f.x), mix( c, d,f.x), f.y);
  
  return 12.0*pow( r, 1.0 );

}

float terrain( vec2 x )
{
  float f = terrainSoft( x );

  float h = smoothstep( 0.4, 0.8, noise( 2.0*x ) );
  f -= 0.2*h;

  float d = noise( 35.0*x.yx*vec2(0.1,1.0) );
  f += 0.003*d * h*h;

  return f;
}


vec2 sdSegment2( vec3 a, vec3 b, vec3 p, float ll )
{
  vec3 pa = p - a;
  vec3 ba = b - a;
  float h = clamp( dot(pa,ba)*ll, 0.0, 1.0 );
  
  return vec2( length( pa - ba*h ), h );
}

vec3 solve( vec3 p, float l1, float l2, vec3 dir )
{
  vec3 q = p*( 0.5 + 0.5*(l1*l1-l2*l2)/dot(p,p) );
  
  float s = l1*l1 - dot(q,q);
  s = max( s, 0.0 );
  q += sqrt(s)*normalize(cross(p,dir));
  
  return q;

}

vec3 solve( vec3 a, vec3 b, float l1, float l2, vec3 dir )
{
  return a + solve( b-a, l1, l2, dir );
}

float smin( float a, float b )
{
#if 0
  float k = 32.0;
  float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
#else
  float k = 0.1;
  float h = clamp( 0.5 + 0.5*(b-a)/k, 0.0, 1.0 );
  return mix( b, a, h ) - k*h*(1.0-h);
#endif  
}

struct Monster
{
  vec3 center;
  vec3 mww;
  vec3 ne[6];
  vec3 f0b[6];
};

  
Monster monster;

vec2 sdMonster( in vec3 p )
{
  vec3 q = p - monster.center;
  
  vec3 muu = vec3(1.0,0.0,0.0);
  vec3 mvv = normalize( cross(monster.mww,muu) );

  q = vec3( q.x, dot(mvv,q), dot(monster.mww,q) );

    // body
  float ab = (0.5 + 0.5*cos( 1.0 + 40.0*pow(0.5-0.5*q.z,2.0) ))*(0.5+0.5*q.z);
  float d1 = length( q*vec3(1.5,2.2,1.0) ) - 1.0 - 0.3*ab;
  d1 += 0.03*sin(20.0*q.z)*(0.5+0.5*clamp(2.0*q.y,0.0,1.0));
  float f = 0.5 - 0.5*q.z;
  d1 += f*0.04*sin(40.0*q.y)*sin(40.0*q.x)*sin(40.0*q.z)*clamp(2.0*q.y+0.5,0.0,1.0);
  float ho = 1.0-clamp( 3.0*abs(q.x), 0.0, 1.0 );
  d1 += 0.1*(1.0-sqrt(1.0-ho*ho))*smoothstep( 0.0,0.1,-q.z );
  
  // legs
  for( int i=0; i<6; i++ )
  {
    float s = -sign( float(i)-2.5 );
    float h = mod( float(i), 3.0 )/3.0;
    
    vec3 bas = monster.center + muu*s*0.5 + monster.mww*1.0*(h-0.33) ;

    vec3 n1 = monster.ne[i];
    vec2 hh = sdSegment2( bas, n1, p, 1.0/(1.6*1.6) );
    d1 = smin( d1, hh.x-mix(0.15,0.05,hh.y) + 0.05*sin(6.2831*hh.y) );
    hh = sdSegment2( n1, monster.f0b[i], p, 1.0/(1.2*1.2) );
    d1 = smin( d1, hh.x-mix(0.06,0.02,hh.y) + 0.01*cos(2.0*6.2831*hh.y) );
  }
  
  
  vec2 res = vec2( 0.5*d1, 1.0 );

  // eyes
  q.x = abs(q.x);
  float d3 = length( q - vec3(0.3,0.05,0.9) ) - 0.3;
  if( d3<res.x ) res = vec2( d3, 0.0 );

  return res;
}


vec2 map( in vec3 p )
{
  // monster
  vec2 res = sdMonster( p );

  // terrain
  float d2 = 0.7*(p.y - terrain(p.xz));
  if( d2<res.x ) res=vec2(d2,2.0);
  
  return res;
}

vec3 intersect( in vec3 ro, in vec3 rd )
{
  float maxd = 70.0;
  float precis = 0.001;
  float h=precis*2.0;
  float t = 0.0;
  float d = 0.0;
  float m = 1.0;
  for( int i=0; i<120; i++ )
  {
    if( abs(h)<precis||t>maxd ) continue;//break;
    t += h;
    vec2 res = map( ro+rd*t );
    h = res.x;
    d = res.y;
    m = res.y;
    // Get less precise as distance grows.
    precis = 0.0005 * t;
  }

  if( t>maxd ) m=-1.0;
  return vec3( t, d, m );
}

vec3 calcNormal( in vec3 pos )
{
  vec3 eps = vec3(0.002,0.0,0.0);

  return normalize( vec3(
           map(pos+eps.xyy).x - map(pos-eps.xyy).x,
           map(pos+eps.yxy).x - map(pos-eps.yxy).x,
           map(pos+eps.yyx).x - map(pos-eps.yyx).x ) );
}

float softshadow( in vec3 ro, in vec3 rd, float mint, float k )
{
  float res = 1.0;
  float t = mint;
  float h = 1.0;
  for( int i=0; i<32; i++ )
  {
    h = map(ro + rd*t).x;
    res = min( res, k*h/t );
    t += clamp( h, 0.07, 1.0 );
  }
  return clamp(res,0.0,1.0);
}


vec3 lig = normalize(vec3(-1.0,0.4,0.2));

vec3 path( float t )
{
  vec3 pos = vec3( 0.0 );
  pos.z += t*0.4;
  pos.y = 1.0 + terrainSoft( pos.xz );
  return pos;
}

// boxplorer 3d hackery
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

void main(void)
{
  vec2 q = gl_FragCoord.xy / iResolution.xy;
  vec2 p = -1.0 + 2.0 * q;
  p.x *= iResolution.x/iResolution.y;
  vec2 m = vec2(0.5);
  // Boxplorer uses mouse for viewing itself just fine.
  //if( iMouse.z>0.0 ) m = iMouse.xy/iResolution.xy;
  m = dir.xy;


  //-----------------------------------------------------
  // animate
  //-----------------------------------------------------
  
  float ctime = 15.0 + iGlobalTime;

  // mobe body
  float atime = 2.0*ctime;
  float ac = noise( 0.5*ctime );
  atime += 4.0*ac;
  monster.center = path( atime );
  vec3 centerN = path( atime+2.0 );
  monster.mww = normalize( centerN - monster.center );
  monster.center.y -= 0.25;

  
  // move legs
  for( int i=0; i<6; i++ )
  {
    float s = -sign( float(i)-2.5 );
    float h = mod( float(i), 3.0 )/3.0;

    float z = 0.5*atime + 1.0*h + 0.25*s;
    float iz = floor(z);
    float fz = fract(z);
    float az = clamp((fz-0.66)/0.34,0.0,1.0);
    
    vec3 fo = vec3(s*1.5, 0.7*az*(1.0-az), (iz + az + (h-0.3)*4.0)*0.4*2.0 );
    fo.y += terrain( fo.xz );
    monster.f0b[i] = fo;
    
    vec3 ba = monster.center + vec3(1.0,0.0,0.0)*s*0.5 + monster.mww*1.0*(h-0.33) ;

    monster.ne[i] = solve( ba, fo, 1.6, 1.2, s*vec3(0.0,0.0,-1.0) );
  }

  
  //-----------------------------------------------------
  // camera
  //-----------------------------------------------------
  
  // follow the monster
  float an = 0.0 + 0.1*ctime - 6.28*m.x;
  float cr = 0.3*cos(0.2*ctime);
  vec3 ro = monster.center + vec3(4.0*sin(an),0.2,4.0*cos(an));
  vec3 ta = monster.center;
  ro.y = 0.5 + terrainSoft( ro.xz );
  
  // shake
  ro += 0.04*sin(4.0*ctime*vec3(1.1,1.2,1.3)+vec3(3.0,0.0,1.0) );
  ta += 0.04*sin(4.0*ctime*vec3(1.7,1.5,1.6)+vec3(1.0,2.0,1.0) );
  
  // camera matrix
  vec3 ww = normalize( ta - ro );
  vec3 uu = normalize( cross(ww,vec3(sin(cr),cos(cr),0.0) ) );
  vec3 vv = normalize( cross(uu,ww));
  
  // barrel distortion  
  float r2 = p.x*p.x*0.32 + p.y*p.y;
  p *= (7.0-sqrt(37.5-11.5*r2))/(r2+1.0);
  
  // create view ray
  vec3 rd = normalize( p.x*uu + p.y*vv + 3.0*ww );

  // Use boxplorer stereo camera instead.
  if (!setup_stereo(ro, rd)) {
    gl_FragColor = vec4(0);
    return;
  }

  //-----------------------------------------------------
  // render
  //-----------------------------------------------------
  float nds = clamp(dot(rd,lig),0.0,1.0);
  vec3 bgc = vec3(0.9+0.1*nds,0.95+0.05*nds,1.0)*(0.7 + 0.3*rd.y)*0.98;
  vec3 col = bgc;

  // raymarch
  vec3 tmat = intersect(ro,rd);
  if( tmat.z>-0.5 )
  {
    // geometry
    vec3 pos = ro + tmat.x*rd;
    vec3 nor = calcNormal(pos);
    vec3 ref = reflect( rd, nor );

    // materials
    float mocc = 1.0;
    vec4 mate = vec4(0.0);
    
    // eyes
    if( tmat.z<0.5 )
    {
      mate = 3.0*vec4(0.002,0.002,0.002,10.0);
      mate.xyz *= 0.8 + 0.2*sin(2.0*ref);
      mate.w *= 0.8 + 0.2*sin(20.0*ref.x)*sin(20.0*ref.y)*sin(20.0*ref.z);
      mate.xyz += 0.005*vec3(1.0,0.1,0.0);
    }
    // body
    else if( tmat.z<1.5 )
    {
      // do shading in mosnter space      
      vec3 q = pos - monster.center;
      vec3 muu = vec3(1.0,0.0,0.0);
      vec3 mvv = normalize( cross(monster.mww,muu) );
      q = vec3( q.x, dot(mvv,q), dot(monster.mww,q) );
      vec3 n = vec3( nor.x, dot(mvv,nor), dot(monster.mww,nor) );

      q.x = abs( q.x );
      // base color     
      mate.xyz = 0.3*vec3(1.0,0.7,0.4);     
      mate = mix( mate, 0.3*vec4(1.0,0.9,0.6,1.0), smoothstep(0.0,1.0,-n.y) );
      
            // texture
      mate.xyz *= 0.4+0.6*sqrt(smoothstep(0.0,0.7,texturize(iChannel0,0.5*q,n).xyz));
      
      // stripes
      float ss = smoothstep( 0.5, 0.8, texture2D( iChannel1, 0.5*q.xz*vec2(1.0,0.1) ).x )*smoothstep( 0.0, 0.3, nor.y );
      mate.xyz += 0.15*ss;

      // color adjustment
      mate.w = 0.5;
      mate.xyz *= 0.8;
      mate.xyz *= 0.2+2.3*mate.xzy*vec3(1.0,1.6,1.4);
      
            // occlusion      
      float ho = 1.0-clamp( 5.0*abs(q.x), 0.0, 1.0 );
      mocc *= 1.0-(1.0-sqrt(1.0-ho*ho))*smoothstep( 0.0,0.1,-q.z );
      
      // bump
      vec3 bnor = -1.0 + 2.0*texturize( iChannel2, 3.0*q, n ).xyz;
      bnor.y = abs(bnor.y);
      nor = normalize( nor + (1.0-ss)*0.2*normalize(bnor) );

    }
        // terrain    
    else if( tmat.z<2.5 )
    {
      mate = vec4(1.0, 0.9, 0.5, 0.0);
      float nn =  noise( 2.0*pos.xz );

      mate.xyz = mix( 0.7*mate.xyz, mate.xyz*0.65*vec3(0.8,0.9,1.0), 1.0-smoothstep( 0.4, 0.9, nn ) );
      mate.xyz *= 0.45;

      vec3 ff  = 0.05+0.95*texture2D( iChannel0, 10.0*0.008*pos.xz ).xyz;
      ff *= 0.05+0.95*texture2D( iChannel0, 10.0*0.052*pos.xz ).xyz;
      ff *= 0.05+0.95*texture2D( iChannel0, 10.0*0.403*pos.xz ).xyz;
      
      float aa = mix( 1.0, 0.3, smoothstep( 0.65, 0.8, nn ) );
      mate.xyz *= (1.0-aa) + aa*sqrt(ff)*3.0;
      
      float d = smoothstep( 0.0, 0.5, abs(nn-0.75) );
      mate.xyz *= 0.6+0.4*d;
      d = smoothstep( 0.0, 0.2, abs(nn-0.75) );
      mocc *= 0.7+0.3*d;
      mocc *= 0.03+0.97*pow( clamp( 0.5*length( pos.xz - monster.center.xz ), 0.0, 1.0 ), 2.0 );

      vec3 bnor = -1.0 + 2.0*texture2D( iChannel2, 3.0*pos.xz ).xyz;
      bnor.y = abs(bnor.y);
      nor = normalize( nor + 0.1*normalize(bnor) );
    }


    // lighting
    float occ = (0.5 + 0.5*nor.y)*mocc;
    float amb = 0.5 + 0.5*nor.y;
    float dif = max(dot(nor,lig),0.0);
    float bac = max(dot(nor,-lig),0.0);
    float sha = 0.0; if( dif>0.0 ) sha=softshadow( pos, lig, 0.05, 32.0 );
    float fre = pow( clamp( 1.0 + dot(nor,rd), 0.0, 1.0 ), 2.0 );
    float spe = pow( clamp( dot(lig,reflect(rd,nor)), 0.0, 1.0), 6.0 );
    
    // lights
    vec3 brdf = vec3(0.0);
    brdf += 3.0*dif*vec3(1.50,1.00,0.65)*pow(vec3(sha),vec3(1.0,1.2,1.5));
    brdf += 5.0*amb*vec3(0.12,0.11,0.10)*occ;
    brdf += 3.0*bac*vec3(0.30,0.20,0.15)*occ;
    brdf += 1.0*fre*vec3(0.50,0.50,0.50)*occ*15.0*mate.w*(0.05+0.95*dif*sha);
    brdf += 1.0*spe*vec3(1.0)*6.0*occ*mate.w*dif*sha;

    // surface-light interacion
    col = mate.xyz* brdf;

    // fog    
    col = mix( col, 0.8*bgc, clamp(1.0-1.2*exp(-0.013*tmat.x ),0.0,1.0) );
  }
  else
  {
    // sun    
    vec3 sun = vec3(1.0,0.8,0.5)*pow( nds, 24.0 );
    col += sun;
  }

  // sun scatter
  col += 0.4*vec3(0.2,0.14,0.1)*pow( nds, 7.0 );


  //-----------------------------------------------------
  // postprocessing
  //-----------------------------------------------------
  // gamma
  col = pow( col, vec3(0.45) );

  // desat
  col = mix( col, vec3(dot(col,vec3(0.333))), 0.2 );

  // tint
  col *= vec3( 1.0, 1.0, 1.0*0.9);

  // vigneting
  col *= 0.2 + 0.8*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.1 );

  // fade in  
  col *= smoothstep( 0.0, 2.0, iGlobalTime );

  gl_FragColor = vec4( col, 1.0 );
}
