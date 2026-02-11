// Final stage: optional zoom
varying vec2 iTexCoord;

uniform sampler2D iTexture;
uniform vec3 iZoom;

void main() {
  vec4 color;

  if (iZoom.z != 1.0) {
    // optional zoom.
    vec2 iZoomCenter = vec2(iZoom.x, iZoom.y);
    float iZoomFactor = clamp(iZoom.z, 1.f, 16.f);
    vec2 iZoomOffset = iTexCoord - iZoomCenter;
    vec2 iZoomSampler = iZoomCenter + iZoomOffset / iZoomFactor;
    vec4 cur = texture2D(iTexture, iZoomSampler);
    // show classy full white live pixels when zoomed in.
    color = vec4(int(cur.a));
  } else {
    color = texture2D(iTexture, iTexCoord);
  }

  gl_FragColor = color;
}
