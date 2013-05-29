// Simple copy

varying vec2 texture_coordinate;

uniform sampler2D my_texture;

void main() {
  gl_FragColor = texture2D(my_texture, texture_coordinate);
}
