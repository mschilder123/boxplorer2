varying vec2 iTexCoord;
void main() {
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  iTexCoord = vec2(gl_MultiTexCoord0);
}
