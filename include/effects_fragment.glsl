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

void main() {
  vec4 color;

  // blur?
  if (enable_dof) {
    color = min(texture2D(iBlur0, iTexCoord),
                texture2D(iBlur1, iTexCoord));
  } else {
    color = texture2D(iTexture, iTexCoord);
  }

  // tonemap
  float YD = exposure * (exposure / maxBright + 1.0) / (exposure + 1.0);
  color *= YD;

  // gamma
  float invGamma = 1.0 / gamma;
  gl_FragColor = pow(color, vec4(invGamma, invGamma, invGamma, 1.0));

  // copy Z
  gl_FragDepth = texture2D(iDepth, iTexCoord).x;
}
