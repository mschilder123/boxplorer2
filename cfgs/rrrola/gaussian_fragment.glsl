// 9-tap Gaussian blur shader
// Two stages, horizontal and vertical.
// Run 2nd stage on output of first.

uniform sampler2D iTexture;  // input texture
uniform int iXorY;  // Whether this is X or Y direction pass
uniform float xres, yres;  // iTexture size

varying vec2 iTexCoord;  // from vertex shader

vec2 texRez = vec2(1.0/xres, 1.0/yres);

const float gb_offset[3] = { 0.0, 1.3846153846, 3.2307692308 };
const float gb_weight[3] = { 0.2270270270, 0.3162162162, 0.0702702703 };

#if 0  // gamma? What gamma
#define P(x) pow((x), vec3(2.2))
#define S(x) pow((x), vec3(1/2.2))
#else
#define P(x) (x)
#define S(x) (x)
#endif

vec3 mixDimension(sampler2D text, vec2 coord) {
  vec3 col = P(texture2D(text, coord).rgb) * gb_weight[0];
  vec2 xory = vec2(iXorY, 1. - iXorY) * texRez;

  for (int i = 1; i < 3; ++i) {
    col += P(texture2D(text, coord + xory * gb_offset[i]).rgb) * gb_weight[i];
    col += P(texture2D(text, coord - xory * gb_offset[i]).rgb) * gb_weight[i];
  }
  return col;
}

void main() {
  vec4 a = texture2D(iTexture, iTexCoord);
#if 0  // visual debug
  if (gl_FragCoord.y < yres/2.0) gl_FragColor = a; else
#endif
  gl_FragColor = vec4(S(mixDimension(iTexture, iTexCoord)), a.a);
}
