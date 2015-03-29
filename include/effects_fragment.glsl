// Final stage: blend blur, bloom, tonemap, gamma.
varying vec2 iTexCoord;
uniform sampler2D iTexture;
uniform sampler2D iDepth;
uniform sampler2D iBlur0;
uniform sampler2D iBlur1;
uniform bool enable_dof;

uniform float exposure;
uniform float maxBright;
uniform float gamma;

#if defined(ST_OCULUS)
uniform float xres;
uniform float yres;
uniform float ipd;
#endif

void main() {
  vec4 color;

#if defined ST_OCULUS
  // Attempt chromatic aberration compensation.
  vec2 uv = gl_FragCoord.xy / vec2(xres, yres);
  float halfIpd = (ipd * 10. - 20.) / xres;

  vec2 o;

  if (uv.x < .5) {
    o = vec2(.25 - halfIpd, .5);
  } else {
    o = vec2(.75 + halfIpd, .5);
  }
  
  vec2 q = uv - o;
  vec2 p = vec2(4., 2.) * q;  // [-1~1]

  float l2 = length(p * vec2(.94, 1.015));  // magic value, eyeballed.
  if (l2 > 1.0) {
    // Don't waste time on pixels we can't see.
    gl_FragColor = vec4(.0);
    gl_FragDepth = 0.;
    return;
  }

  // Chromatic aberration pre-compensation, eyeballed.
  // TODO: nonlinear?
  const float Raber = 1.000;
  const float Gaber = 1.013;
  const float Baber = 1.032;
  p *= vec2(.25, .5);
  color = vec4(
    texture2D(iTexture, o + p * Raber).r,
    texture2D(iTexture, o + p * Gaber).g,
    texture2D(iTexture, o + p * Baber).b,
    1.0);

  // Vignette border, outer 10%.
  float dropOff = clamp(l2 - .9, 0., .1) * 10.;
  color *= (1. - dropOff);
 
#else  // ST_OCULUS

  // blur
  if (enable_dof) {
    color = min(texture2D(iBlur0, iTexCoord),
                texture2D(iBlur1, iTexCoord));
  } else {
    color = texture2D(iTexture, iTexCoord);
  }

#endif  // ST_OCULUS

  // tonemap
  float Y = dot(vec3(0.30, 0.59, 0.11), color.rgb);
  // TODO: use Y?
  float YD = exposure * (exposure / maxBright + 1.0) / (exposure + 1.0);
  color *= YD;

  // gamma
  float invGamma = 1.0 / gamma;
  gl_FragColor = pow(color, vec4(invGamma, invGamma, invGamma, 1.0));

  // copy Z
  gl_FragDepth = texture2D(iDepth, iTexCoord).x;
}
