#pragma once

#include <math.h>
#include <string>

#include "interpolate.h"
#include "params.h"
#include "shader_procs.h"
#include "uniforms.h"
#include "utils.h"

#if defined(_WIN32)
#pragma warning(disable : 4996) // unsafe function
#pragma warning(disable : 4244) // double / float conversion
#pragma warning(disable : 4305) // double / float truncation
#pragma warning(disable : 4800) // forcing value to bool
#endif                          // _WIN32

extern std::string de_func_name;
extern Uniforms uniforms;

// Pinhole camera modes.
enum StereoMode {
  ST_NONE = 0,
  ST_OVERUNDER,
  ST_XEYED,
  ST_INTERLACED,
  ST_SIDEBYSIDE,
  ST_QUADBUFFER,
  ST_OCULUS,
  ST_ANAGLYPH,
  ST_SPHERICAL,
  ST_DOME,
  ST_COMPUTE_DE_ONLY
};

// For seamless cube rendering.
typedef int ViewQuadrant;
#define VQ_FRONT 0
#define VQ_BACK 1
#define VQ_UP 2
#define VQ_DOWN 3
#define VQ_RIGHT 4
#define VQ_LEFT 5
#define VQ_DONE 6
#define VQ_LETTER "fbudrl"

class KeyFrame {
public:
  double v[16]; // view matrix
  double q[4];  // quaterion orientation
  double x[4];  // r4 splineable q

  // Declare common parameters.
#define PROCESS(a, b, c, d) a b;
  PROCESS_COMMON_PARAMS
#undef PROCESS

  // Par[] parameter array.
  float par[NUMPARS][3]; // min(|this|, |glsl|) gets sent to shader.

  // Shader other uniforms.
  int iunis[NUMPARS];
  float funis[NUMPARS];

  int n_iunis;
  int n_funis;
  bool isKey_; // Whether this frame is actually a defined KeyFrame.

  KeyFrame();

  double distanceTo(const KeyFrame &other) const;

  double *right() { return &v[0]; }
  double *up() { return &v[4]; }
  double *ahead() { return &v[8]; }
  double *pos() { return &v[12]; }

  void setKey(bool key) { isKey_ = key; }
  bool isKey() const { return isKey_; }

  // Orthogonalize v[]
  void orthogonalize();

  // Rotate q[] by `deg` degrees around a normalized axis
  // and set rotation part of v[] to q[].
  // Behaves like `glRotate` without normalizing the axis.
  void rotate(double deg, double x, double y, double z);

  // Move camera in a direction relative to the view direction.
  // Behaves like `glTranslate`.
  void move(double x, double y, double z) {
    for (int i = 0; i < 3; i++) {
      pos()[i] += right()[i] * x + up()[i] * y + ahead()[i] * z;
    }
    iBackbufferCount = 0;
  }

  // Move camera in the normalized absolute direction `dir` by `len` units.
  void moveAbsolute(double *dir, double len) {
    for (int i = 0; i < 3; i++) {
      pos()[i] += len * dir[i];
    }
  }

  void set(double x, double y, double z) {
    pos()[0] = x;
    pos()[1] = y;
    pos()[2] = z;
  }

  // Map a uniform name to a address within this.
  // Returns NULL on fail.
  void *map_address(const std::string &type, const std::string &name, int n);
};

class Camera : public KeyFrame {
public:
  Camera &operator=(const KeyFrame &other) {
    *((KeyFrame *)this) = other;
    iBackbufferCount = 0; // reset progressive rendering count.
    return *this;
  }

  // Set the OpenGL modelview matrix to the camera matrix, for shader.
  void activate(ViewQuadrant vq) {
    orthogonalize();
    glMatrixMode(GL_MODELVIEW);
    // Tweak view for desired quadrant.
    double m[16];
    memcpy(m, v, sizeof(m));
    switch (vq) {
    case VQ_FRONT:
      // m[0] = v[0]; m[1] = v[1]; m[2] = v[2];  // right
      // m[4] = v[4]; m[5] = v[5]; m[6] = v[6];  // up
      // m[8] = v[8]; m[9] = v[9]; m[10] = v[10];  // ahead
      break;
    case VQ_BACK:
      // ahead = -ahead; right = -right
      m[0] = -v[0];
      m[1] = -v[1];
      m[2] = -v[2];
      m[8] = -v[8];
      m[9] = -v[9];
      m[10] = -v[10];
      break;
    case VQ_UP:
      // ahead = up; up = -ahead;
      m[4] = -v[8];
      m[5] = -v[9];
      m[6] = -v[10];
      m[8] = v[4];
      m[9] = v[5];
      m[10] = v[6];
      break;
    case VQ_DOWN:
      // ahead = -up; up = ahead
      m[4] = v[8];
      m[5] = v[9];
      m[6] = v[10];
      m[8] = -v[4];
      m[9] = -v[5];
      m[10] = -v[6];
      break;
    case VQ_RIGHT:
      // ahead = right; right = -ahead
      m[0] = -v[8];
      m[1] = -v[9];
      m[2] = -v[10];
      m[8] = v[0];
      m[9] = v[1];
      m[10] = v[2];
      break;
    case VQ_LEFT:
      // ahead = -right; right = ahead
      m[0] = v[8];
      m[1] = v[9];
      m[2] = v[10];
      m[8] = -v[0];
      m[9] = -v[1];
      m[10] = -v[2];
      break;
    }
    glLoadMatrixd(m);
  }

  // Set the OpenGL modelview and projection for gl*() functions.
  void activateGl() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double z_near = fabs(speed);
    double z_far = speed * 65535.0;
    double fH = tan(fov_y * PI / 360.0f) * z_near;
    double fW = tan(fov_x * PI / 360.0f) * z_near;
    glFrustum(-fW, fW, -fH, fH, z_near, z_far);

    orthogonalize();
    double matrix[16] = {right()[0], up()[0], -ahead()[0], 0,
                         right()[1], up()[1], -ahead()[1], 0,
                         right()[2], up()[2], -ahead()[2], 0,
                         0,          0,       0,           1};
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(matrix);
    // Do not translate, keep eye at 0 to retain drawing precision.
    // glTranslated(-pos()[0], -pos()[1], -pos()[2]);
  }

  // Load configuration.
  bool loadConfig(const std::string &configFile, std::string *defines = NULL) {
    bool result = false;
    std::string filename(WorkingDir + configFile);
    FILE *f;
    if ((f = fopen(filename.c_str(), "r")) != 0) {
      size_t i;
      char s[32768];                     // max line length
      while (fscanf(f, " %s", s) == 1) { // read word
        if (s[0] == 0 || s[0] == '#')
          continue;

        int v;

        // Parse #defines out of config.cfg to prepend to .glsl
        if (defines) {
          if (!strcmp(s, "d") || !strcmp(s, "c")) {
            std::string a(s);
            v = fscanf(f, " %s", s);
            if (v == 1) {
              std::string define = "#define " + a + " " + s + "\n";
              printf("%s : %s", __func__, define.c_str());
              defines->append(define);
              if (!a.compare("d")) {
                de_func_name.assign(s);
                printf("%s : de_func %s\n", __func__, de_func_name.c_str());
              }
            }
          }
        }

        double val;

        if (!strcmp(s, "position")) {
          v = fscanf(f, " %lf %lf %lf", &pos()[0], &pos()[1], &pos()[2]);
          continue;
        }
        if (!strcmp(s, "direction")) {
          v = fscanf(f, " %lf %lf %lf", &ahead()[0], &ahead()[1], &ahead()[2]);
          continue;
        }
        if (!strcmp(s, "upDirection")) {
          v = fscanf(f, " %lf %lf %lf", &up()[0], &up()[1], &up()[2]);
          continue;
        }

        // Parse common parameters.
#define PROCESS(type, name, nameString, doSpline)                              \
  if (!strcmp(s, nameString)) {                                                \
    v = fscanf(f, " %lf", &val);                                               \
    name = val;                                                                \
    continue;                                                                  \
  }
        PROCESS_COMMON_PARAMS
#undef PROCESS

        for (i = 0; i < ARRAYSIZE(par); i++) {
          char p[256];
          sprintf(p, "par%lu", (unsigned long)i);
          if (!strcmp(s, p)) {
            v = fscanf(f, " %f %f %f", &par[i][0], &par[i][1], &par[i][2]);
            break;
          }
        }
      }
      fclose(f);
      printf("%s : read '%s'\n", __func__, configFile.c_str());
      result = true;
    } else {
      printf("%s : failed to open '%s'\n", __func__, configFile.c_str());
    }
    if (result)
      sanitizeParameters();
    return result;
  }

  // Make sure parameters are OK.
  void sanitizeParameters(void) {
    // Resolution: if only one coordinate is set, keep 4:3 aspect ratio.
    if (width < 1) {
      if (height < 1) {
        height = 480;
      }
      width = height * 4 / 3;
    }
    if (height < 1)
      height = width * 3 / 4;

    // FOV: keep pixels square unless stated otherwise.
    // Default FOV_y is 42 degrees.
    if (fov_x <= 0) {
      if (fov_y <= 0) {
        // ~60x42 degrees is fov for normal position in front of a monitor.
        fov_y = 42;
      }
      fov_x = atan(tan(fov_y * PI / 180 / 2) * width / height) / PI * 180 * 2;
    }
    if (fov_y <= 0)
      fov_y = atan(tan(fov_x * PI / 180 / 2) * height / width) / PI * 180 * 2;

    // Fullscreen: 0=off, anything else=on.
    if (fullscreen != 0 && fullscreen != 1)
      fullscreen = 1;

    // The others are easy.
    if (multisamples < 1)
      multisamples = 1;
    // if (speed <= 0) speed = 0.005;  // units/frame
    if (keyb_rot_speed <= 0)
      keyb_rot_speed = 5; // degrees/frame
    if (mouse_rot_speed <= 0)
      mouse_rot_speed = 1; // degrees/pixel
    if (max_steps < 1)
      max_steps = 128;
    if (min_dist <= 0)
      min_dist = 0.0001;
    if (iters < 1)
      iters = 13;
    if (color_iters < 0)
      color_iters = 9;
    if (ao_eps <= 0)
      ao_eps = 0.0005;
    if (ao_strength <= 0)
      ao_strength = 0.1;
    if (glow_strength <= 0)
      glow_strength = 0.25;
    if (dist_to_color <= 0)
      dist_to_color = 0.2;

    if (exposure == 0)
      exposure = 1.0;
    if (maxBright == 0)
      maxBright = 1.0;
    if (gamma == 0)
      gamma = 1.0;

    iBackbufferCount = 0; // No samples in backbuffer yet.

    orthogonalize();
    mat2quat(this->v, this->q);

    // Don't do anything with user parameters - they must be
    // sanitized (clamped, ...) in the shader.
  }

  // Save configuration.
  void saveConfig(const std::string &configFile, std::string *defines = NULL) {
    FILE *f;
    std::string filename(WorkingDir + configFile);
    if ((f = fopen(filename.c_str(), "w")) != 0) {
      if (defines != NULL)
        fprintf(f, "%s", defines->c_str());

        // Write common parameters.
#define PROCESS(type, name, nameString, doSpline)                              \
  fprintf(f, nameString " %g\n", (double)name);
      PROCESS_COMMON_PARAMS
#undef PROCESS

      fprintf(f, "position %12.12e %12.12e %12.12e\n", pos()[0], pos()[1],
              pos()[2]);
      fprintf(f, "direction %g %g %g\n", ahead()[0], ahead()[1], ahead()[2]);
      fprintf(f, "upDirection %g %g %g\n", up()[0], up()[1], up()[2]);
      for (size_t i = 0; i < ARRAYSIZE(par); i++) {
        fprintf(f, "par%lu %g %g %g\n", (unsigned long)i, par[i][0], par[i][1],
                par[i][2]);
      }
      fclose(f);
      printf("%s : wrote '%s'\n", __func__, filename.c_str());
    }
  }

  // Send parameters to gpu.
  void setUniforms(float x_scale, float x_offset, float y_scale, float y_offset,
                   double spd, GLuint program) {
#define glSetUniformf(name)                                                    \
  glUniform1f(glGetUniformLocation(program, #name), name);
#define glSetUniformfv(name)                                                   \
  glUniform3fv(glGetUniformLocation(program, #name), ARRAYSIZE(name),          \
               (float *)name);
#define glSetUniformi(name)                                                    \
  glUniform1i(glGetUniformLocation(program, #name), name);

    // These might be dupes w/ uniforms.send() below.
    // Leave for now until all .cfg got updated.
    glSetUniformi(max_steps);
    glSetUniformf(min_dist);
    glSetUniformi(iters);
    glSetUniformi(color_iters);
    glSetUniformf(ao_eps);
    glSetUniformf(ao_strength);
    glSetUniformf(glow_strength);
    glSetUniformf(dist_to_color);
    glSetUniformi(nrays);
    glSetUniformf(focus);

    // Non-user uniforms.
    glSetUniformf(fov_x);
    glSetUniformf(fov_y);

    glSetUniformf(x_scale);
    glSetUniformf(x_offset);
    glSetUniformf(y_scale);
    glSetUniformf(y_offset);

    glSetUniformf(time);

    glUniform1f(glGetUniformLocation(program, "speed"), spd);
    glUniform1f(glGetUniformLocation(program, "ipd"), ipd);
    glUniform1f(glGetUniformLocation(program, "xres"), width);
    glUniform1f(glGetUniformLocation(program, "yres"), height);

    // Also pass in some double precision values, if supported.
    if (glUniform1d) {
      glUniform1d(glGetUniformLocation(program, "dspeed"), spd);
      // For some reason 3dv below stopped working reliably..
      glUniform1d(glGetUniformLocation(program, "deyex"), pos()[0]);
      glUniform1d(glGetUniformLocation(program, "deyey"), pos()[1]);
      glUniform1d(glGetUniformLocation(program, "deyez"), pos()[2]);
    }
    if (glUniform3dv) {
      glUniform3dv(glGetUniformLocation(program, "deye"), 3, pos());
    }

    // Old-style par[] list.
    glSetUniformfv(par);

#undef glSetUniformf
#undef glSetUniformfv
#undef glUniform1i

    // New-style discovered && active uniforms only.
    uniforms.send(program);
  }

  void render(enum StereoMode stereo, ViewQuadrant vq, GLuint program,
              int polarity) {
    activate(vq); // Load view matrix for shader.
    switch (stereo) {
    case ST_OVERUNDER: { // left / right
      setUniforms(1.0, 0.0, 2.0, 1.0, +speed, program);
      glRects(-1, -1, 1, 0); // draw bottom half of screen
      setUniforms(1.0, 0.0, 2.0, -1.0, -speed, program);
      glRects(-1, 0, 1, 1); // draw top half of screen
    } break;
    case ST_QUADBUFFER: { // left - right
      glDrawBuffer(GL_BACK_LEFT);
      setUniforms(1.0, 0.0, 1.0, 0.0, -speed * polarity, program);
      glRects(-1, -1, 1, 1);
      glDrawBuffer(GL_BACK_RIGHT);
      setUniforms(1.0, 0.0, 1.0, 0.0, +speed * polarity, program);
      glRects(-1, -1, 1, 1);
    } break;
    case ST_XEYED: { // right | left
      setUniforms(2.0, +1.0, 1.0, 0.0, +speed, program);
      glRectf(-1, -1, 0, 1); // draw left half of screen
      setUniforms(2.0, -1.0, 1.0, 0.0, -speed, program);
      glRectf(0, -1, 1, 1); // draw right half of screen
    } break;
    case ST_SIDEBYSIDE: { // left | right
      setUniforms(2.0, +1.0, 1.0, 0.0, -speed, program);
      glRectf(-1, -1, 0, 1); // draw left half of screen
      setUniforms(2.0, -1.0, 1.0, 0.0, +speed, program);
      glRectf(0, -1, 1, 1); // draw right half of screen
    } break;
    case ST_NONE:
    case ST_SPHERICAL:
    case ST_DOME:
      setUniforms(1.0, 0.0, 1.0, 0.0, speed, program);
      glRects(-1, -1, 0, 1); // draw left half
      glRects(0, -1, 1, 1);  // draw right half
      break;
    case ST_INTERLACED:
    case ST_ANAGLYPH:
      setUniforms(1.0, 0.0, 1.0, 0.0, speed * polarity, program);
      glRects(-1, -1, 0, 1); // draw left half
      glRects(0, -1, 1, 1);  // draw right half
      break;
    case ST_OCULUS:
      setUniforms(2.0, +1.0, 1.0, 0.0, -speed, program);
      glRectf(-1, -1, 0, 1); // draw left half of screen
      setUniforms(2.0, -1.0, 1.0, 0.0, +speed, program);
      glRectf(0, -1, 1, 1); // draw right half of screen
      break;
    case ST_COMPUTE_DE_ONLY:
      setUniforms(1.0, 0.0, 1.0, 0.0, speed, program);
      float xrange = 4.0 / width; // Aim for ~8x8 pixels.
      float yrange = 4.0 / height;
      glRectf(-1, 1, -1 + xrange, 1 - yrange); // Only care about top corner.
      break;
    }
  }

  void mixHydraOrientation(float *quat) {
    double q[4];
    q[0] = quat[0];
    q[1] = quat[1];
    q[2] = quat[2];
    q[3] = quat[3];
    qnormalize(q);
    qmul(q, this->q);
    quat2mat(q, this->v);
  }

  // take this->q and q and produce this->v,q := this->q + q
  void mixSensorOrientation(float q[4]) {
    double q1[4];
    q1[0] = q[0];
    q1[1] = q[1];
    q1[2] = q[2];
    q1[3] = q[3];

    q1[2] = -q1[2]; // We roll other way
    qnormalize(q1);

    // combine current view quat with sensor quat.
    qmul(q1, this->q);

    quat2mat(q1, this->v);

    // this->q = q1;
    this->q[0] = q1[0];
    this->q[1] = q1[1];
    this->q[2] = q1[2];
    this->q[3] = q1[3];
  }

  // take this->v and q and produce this->v,q := this->v - q
  void unmixSensorOrientation(float q[4]) {
    double q1[4];
    q1[0] = q[0];
    q1[1] = q[1];
    q1[2] = q[2];
    q1[3] = q[3];

    q1[2] = -q1[2]; // We roll other way
    qnormalize(q1);

    // apply inverse current view.
    qinvert(q1, q1);

    mat2quat(this->v, this->q);
    qmul(q1, this->q);

    quat2mat(q1, this->v);

    // this->q = q1;
    this->q[0] = q1[0];
    this->q[1] = q1[1];
    this->q[2] = q1[2];
    this->q[3] = q1[3];
  }
};
