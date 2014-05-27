// http://ivizlab.sfu.ca/media/DiPaolaMcIntoshRiecke2012.pdf
// Two stages, horizontal and vertical.
// Run 2nd stage on output of first.

uniform sampler2D iTexture;  // input texture
uniform sampler2D iDepth;  // input depth
uniform int iBlurLevel;
uniform int iXorY;  // Whether this is X or Y direction pass
uniform float xres, yres;  // iTexture size

varying vec2 iTexCoord;  // from vertex shader

vec2 texRez = vec2(1.0/xres, 1.0/yres);

// box
//const float gb_offset[13] = float[](-11.5, -9.5, -7.5, -5.5, -3.5, -1.5, 0.0, 1.5, 3.5, 5.5, 7.5, 9.5, 11.5);

#define SQRT2 0.70710678118  // sqrt(2) * .5


#if 1  // gamma? What gamma
#define GAMMA 4.2
#define P(x) pow((x), vec4(GAMMA,GAMMA,GAMMA,1.0))
#define S(x) pow((x), vec4(1.0/GAMMA,1.0/GAMMA,1.0/GAMMA,1.0))
#else
#define P(x) (x)
#define S(x) (x)
#endif

vec4 mixDimension(vec2 coord) {
  const float bleedingBias = 0.02;
  const float bleedingMult = 30.0;

  float totalWeight = 0.0;
  vec4 color = vec4(0.0, 0.0, 0.0, 0.0);

  vec4 centerPixel = P(texture2D(iTexture, coord));
  float centerDepth = texture2D(iDepth, coord).x;
 
  vec2 xory;

  if (iBlurLevel == 0) {
    // (1, 0) or (0, 1) : square
    xory = vec2(float(iXorY), 1. - float(iXorY)) * texRez * centerPixel.a;
  } else {
    // (s2, s2) or (s2, -s2) : square rotate 45 degrees
    xory = vec2(SQRT2, 2.*(.5 - float(iXorY)) * SQRT2) * texRez * centerPixel.a;
  }

  float gb_offset = -11.5;
  for (int i = 0; i < 13; ++i) {
    vec2 sampleCoord = coord + xory * gb_offset;
    gb_offset += 2.0;
    vec4 samplePixel = P(texture2D(iTexture, sampleCoord));
    float sampleDepth = texture2D(iDepth, sampleCoord).x;

    float weight = (sampleDepth < centerDepth) ?
        samplePixel.a * bleedingMult : 1.0;
    weight = (centerPixel.a > samplePixel.a + bleedingBias) ?
        weight : 1.0;
    weight = clamp(weight, 0.0, 1.0);

    color += samplePixel * weight;
    totalWeight += weight;
  }

  color /= totalWeight;
  //return vec4(color.rgb, centerPixel.a);
  return color;
}

void main() {
  vec4 col;

  col = S(mixDimension(iTexCoord));

  gl_FragColor = col;
}
