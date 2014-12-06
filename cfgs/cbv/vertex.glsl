varying vec3 dir;

void main() {
  gl_Position = gl_Vertex;
  dir = vec3(gl_ModelViewMatrix * vec4(0, 0, 1, 0));
}
