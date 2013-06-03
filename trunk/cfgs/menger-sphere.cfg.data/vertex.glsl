varying vec3 eye, dir;

uniform float fov_x, fov_y;  // Field of vision.

// stereoscopy params:
// x_scale, x_offset, y_scale, y_offset
// 1,        0,       1,       0      : normal view
// 1,        0,       2,       1      : bottom half, 50% vertical pixels
// 1,        0,       2,      -1      : top half, 50% vertical pixels
// 2,        1,       1,       0      : left half, 50% x fov
// 2,       -1,       1,       0      : right half, 50% x fov
uniform float x_offset, x_scale, y_offset, y_scale;

float fov2scale(float fov) { return tan(radians(fov/2.0)); }

// Draw an untransformed rectangle covering the whole screen.
// Get camera position and interpolated directions from the modelview matrix.
void main() {
  gl_Position = gl_Vertex;
  eye = vec3(gl_ModelViewMatrix[3]);
  dir = vec3(gl_ModelViewMatrix *
    vec4(
      fov2scale(fov_x)*(x_scale*gl_Vertex.x + x_offset),
      fov2scale(fov_y)*(y_scale*gl_Vertex.y + y_offset),
      1,
      0));
}
