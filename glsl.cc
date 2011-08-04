/*
 * Little test program to experiment w/ compiling glsl code as C++.
 */
#include <stdio.h>
#include <math.h>
#include <time.h>

#if !defined(_WIN32)
#define __FUNCTION__ "glsl"
#else

#pragma warning(disable: 4996) // unsafe function
#pragma warning(disable: 4244) // conversion loss
#pragma warning(disable: 4305) // truncation
#pragma warning(disable: 4800) // forcing value to bool
#endif

#include "glsl.h"
#include "TGA.h"

namespace GLSL {

#include "cfgs/menger.cfg.data/fragment.glsl"

int main(int argc, char* argv[]) {
  printf("glsl as C++ test:\n");

  // Setup up params.
  // TODO: read from .cfg
  iters = 6;
  min_dist = 0.0001;
  max_steps = 1162;
  ao_eps = 0.000792447;
  ao_strength = 0.0707946;
  glow_strength = 5.5;
  par[0] = vec3(0.25, -1.77, 1.59);
  par[1] = vec3(6, 1, 0);
  par[8] = vec3(3.63, 2.64, 1.71);
  par[9] = vec3(0, 0, 0);
  par[11] = vec3(0.95, 0.64, 0.1);
  par[14] = vec3(1.0, 0.8, 0.4);
  par[15] = vec3(.03, .4, .4);
  par[16] = vec3(0, 0, 0);

  float fov_x = 91.3085;
  float fov_y = 75;

  // eye position
  vec3 pos(-0.885181,0.633865,-0.915299);
  // view direction
  vec3 ahead(0.174267,-0.820984,-0.543707);
  vec3 up(-0.976163,-0.0714907,-0.204927);
  vec3 right = up.cross(ahead);

  time_t start, end;
  float minTotal = MAX_DIST;
  float maxTotal = 0;
  int maxSteps = 0;

#define XRES 720
#define YRES 480

  TGA tga(XRES, YRES);

  mat4 proj(right.x, up.x, ahead.x, 0,
            right.y, up.y, ahead.y, 0,
            right.z, up.z, ahead.z, 0,
            0,       0,    0,       1);

  start = clock();
  for (int scr_y = 0; scr_y < YRES; ++scr_y) {
    for (int scr_x = 0; scr_x < XRES; ++scr_x) {
      float dx = -1 + scr_x * (2.0/XRES);  // -1..1
      float dy = -1 + scr_y * (2.0/YRES);  // -1..1

      dx *= tan(radians(fov_x * .5));
      dy *= tan(radians(fov_y * .5));

      vec4 ddir = proj * vec4(dx,dy,1,0);
      vec3 dir(ddir.x, ddir.y, ddir.z);

      float totalD = d(pos);
      float side = sign(totalD);
      totalD *= side * .5;

      vec3 dp = normalize(dir);

      float m_zoom = tan(radians(fov_x * .5)) * length(dir) * .5 / XRES;
      float m_dist = max(min_dist, m_zoom * totalD);
      int steps = rayMarch(pos, dp, totalD, side, m_dist, m_zoom);
      if (totalD < MAX_DIST) {
        vec3 p = pos + dp*totalD;
        vec3 n = normal(p, m_dist * .5);
        vec3 col = rayColor(p,dp,n,totalD,m_dist);
        col = mix(col, glowColor, float(steps)/float(max_steps)*glow_strength);
        tga.set(scr_x,scr_y,col);
      }

      if (totalD < minTotal) minTotal = totalD;
      if (totalD > maxTotal) maxTotal = totalD;
      if (steps > maxSteps) maxSteps = steps;
    }
  }
  end = clock();

  printf("%lf sec\n", (double)(end-start)/CLOCKS_PER_SEC);
  printf("min dist %f\n", minTotal);
  printf("max dist %f\n", maxTotal);
  printf("max steps %d\n", maxSteps);

  if (tga.writeFile("test.tga"))
    printf("wrote ./test.tga\n");

  return 0;
}

} // namespace GLSL

int main(int argc, char* argv[]) {
  return GLSL::main(argc, argv);
}
