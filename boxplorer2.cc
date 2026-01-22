// Standard C/C++
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

// debug
#define XXX_SDL_SetCursor(x)                                                   \
  do {                                                                         \
    fprintf(stderr, "line %d: %s\n", __LINE__, #x);                            \
    SDL_SetCursor(x);                                                          \
  } while (0)

#if !defined(_WIN32)

#include <unistd.h>
#define _strdup strdup
#define __FUNCTION__ "boxplorer2"
#define MAX_PATH 256

#else // _WIN32

#pragma warning(disable : 4996) // unsafe function
#pragma warning(disable : 4244) // double / float conversion
#pragma warning(disable : 4305) // double / float truncation
#pragma warning(disable : 4800) // forcing value to bool

#pragma comment(lib, "SDL2.lib")
#pragma comment(lib, "SDL2main.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "shell32.lib")
#pragma comment(lib, "comdlg32.lib")

#include "oculus_sdk4.h"

#if defined(HYDRA)
#include <sixense.h>
#pragma comment(lib, "sixense.lib")
#pragma comment(lib, "sixense_utils.lib")
#endif // HYDRA

#endif // _WIN32

using namespace std;

// External Libraries
#define NO_SDL_GLEXT
#include <AntTweakBar.h>
#include <SDL.h>
#include <SDL_main.h>
#include <SDL_mouse.h>
#include <SDL_opengl.h>
#include <SDL_thread.h>

// Local Headers
#include "TGA.h"
#include "camera.h"
#include "default_shaders.h"
#include "glsl.h"
#include "interpolate.h"
#include "params.h"
#include "shader.h"
#include "shader_manager.h"
#include "shader_procs.h"
#include "uniforms.h"
#include "utils.h"

#define DEFAULT_CONFIG_FILE "default.cfg"
#define DEFAULT_CONFIG "cfgs/rrrola/" DEFAULT_CONFIG_FILE
#define VERTEX_SHADER_FILE "vertex.glsl"
#define FRAGMENT_SHADER_FILE "fragment.glsl"
#define EFFECTS_VERTEX_SHADER_FILE "effects_vertex.glsl"
#define EFFECTS_FRAGMENT_SHADER_FILE "effects_fragment.glsl"
#define DOF_VERTEX_SHADER_FILE "dof_vertex.glsl"
#define DOF_FRAGMENT_SHADER_FILE "dof_fragment.glsl"
#define FXAA_VERTEX_SHADER_FILE "fxaa_vertex.glsl"
#define FXAA_FRAGMENT_SHADER_FILE "fxaa_fragment.glsl"

// Hackery to get the list of DE and COLORING funcs from the glsl.
map<string, float (*)(GLSL::vec3)> DE_funcs;
map<string, double (*)(GLSL::dvec3)> DE64_funcs;
map<string, GLSL::vec3 (*)(GLSL::vec3)> COLOR_funcs;

class DE_initializer {
public:
  DE_initializer(string name, float (*func)(GLSL::vec3)) {
    printf("DECLARE_DE(%s)\n", name.c_str());
    DE_funcs[name] = func;
  }
  DE_initializer(string name, double (*func)(GLSL::dvec3)) {
    printf("DECLARE_DE(%s)\n", name.c_str());
    // Strip _64 from name.
    size_t x64 = name.find("_64");
    if (x64 != string::npos)
      name.erase(x64);
    DE64_funcs[name] = func;
  }
};
#define DECLARE_DE(a) DE_initializer _init##a(#a, &a);
#define DECLARE_COLORING(a) // not interested in color funcs here

string de_func_name;
double (*de_func_64)(GLSL::dvec3) = NULL;
float (*de_func)(GLSL::vec3) = NULL;

namespace GLSL {

// 'globals' capturing the fragment shader output or providing context.
float gl_FragDepth;
vec4 gl_FragColor;
vec4 gl_FragCoord;

// In the c++ version, these are func ptrs, not straight #defines.
// We assign them based on values in .cfg
float (*d)(vec3);
vec3 (*c)(vec3);

// Compile the fragment shader right here.
// This defines a bunch more 'globals' and functions.
#define ST_NONE
#include "cfgs/menger.cfg.data/fragment.glsl"
#undef ST_NONE

} // namespace GLSL

#define FPS_FRAMES_TO_AVERAGE 20

const char *kKEYFRAME = "keyframe";

void clearGlContext(); // fwd decl.

// Our OpenGL SDL2 window.
static const int kMAXDISPLAYS = 6;
class GFX {
public:
  GFX()
      : display_(-1), width_(0), height_(0), last_x_(SDL_WINDOWPOS_CENTERED),
        last_y_(SDL_WINDOWPOS_CENTERED), last_width_(0), last_height_(0),
        fullscreen_(false), window_(NULL), glcontext_(0) {}

  void init() {
    // requires SDL_Init(SDL_INIT_VIDEO);
    for (ndisplays_ = 0;
         ndisplays_ < SDL_GetNumVideoDisplays() && ndisplays_ < kMAXDISPLAYS;
         ++ndisplays_) {
      SDL_GetCurrentDisplayMode(ndisplays_, &mode_[ndisplays_]);
      SDL_GetDisplayBounds(ndisplays_, &rect_[ndisplays_]);
    }
  }

  ~GFX() { reset(); }

  void reset() {
    if (glcontext_) {
      clearGlContext();
      SDL_GL_DeleteContext(glcontext_);
      glcontext_ = 0;
    }
    if (window_) {
      SDL_DestroyWindow(window_);
      window_ = NULL;
    }
    display_ = -1;
  }

  bool resize(int w, int h) {
    int d = display_;

    // ignore resize events when fullscreen.
    if (fullscreen_)
      return false;

    // ignore resize events for fullscreen width.
    if (d != -1 && w == rect_[d].w)
      return false;

    if (window_) {
      // capture current display.
      d = SDL_GetWindowDisplayIndex(window_);
      // capture current window position.
      SDL_GetWindowPosition(window_, &last_x_, &last_y_);
    }

    reset();

    printf(__FUNCTION__ ": %dx%d display %d\n", w, h, d);
    window_ = SDL_CreateWindow("test", last_x_, last_y_, w, h,
                               SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    glcontext_ = SDL_GL_CreateContext(window_);
    display_ = d;
    last_width_ = width_ = w;
    last_height_ = height_ = h;
    return enableShaderProcs(); // re-fetch ptrs in new context.
  }

  bool toggleFullscreen() {
    if (!window_)
      return false;

    // capture current display.
    int d = SDL_GetWindowDisplayIndex(window_);

    if (!fullscreen_) {
      // capture current window position.
      SDL_GetWindowPosition(window_, &last_x_, &last_y_);
    }

    reset();

    // Oculus / Acer hackery:
    // if current resolution matches a display, go
    // fullscreen on that display.
    // Otherwise, stick with current.
    bool likelyOculus = (height_ == 800);
    bool foundMatch = false;
    int alternate1080p = -1;
    for (int i = 0; i < ndisplays_; ++i) {
      if ((width_ == rect_[i].w && height_ == rect_[i].h)) {
        d = i; // exact match, done.
        foundMatch = true;
        break;
      }
      if (rect_[i].h == 1080)
        alternate1080p = i;
      printf("screen %d: %dx%d\n", i, rect_[i].w, rect_[i].h);
    }

    int targetWidth = rect_[d].w;
    int targetHeight = rect_[d].h;

    if (!foundMatch && likelyOculus && alternate1080p != -1) {
      // Could not find exact match.
      // Oculus might be duplicating a 1080p desktop.
      d = alternate1080p;
      // Stick w/ oculus resolution, rather than native screen one.
      targetWidth = width_;
      targetHeight = height_;
    }

    if (!fullscreen_) {
      printf(__FUNCTION__ ": to fullscreen %dx%d display %d\n", targetWidth,
             targetHeight, d);
      window_ = SDL_CreateWindow(
          "boxplorer2", rect_[d].x, rect_[d].y, targetWidth, targetHeight,
          SDL_WINDOW_OPENGL | SDL_WINDOW_FULLSCREEN_DESKTOP);
      width_ = targetWidth;
      height_ = targetHeight;
    } else {
      printf(__FUNCTION__ ": from fullscreen %dx%d display %d\n", last_width_,
             last_height_, d);
      window_ = SDL_CreateWindow("boxplorer2", last_x_, last_y_, last_width_,
                                 last_height_,
                                 SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
      width_ = last_width_;
      height_ = last_height_;
    }

    glcontext_ = SDL_GL_CreateContext(window_);
    display_ = d;

    fullscreen_ = !fullscreen_;
    return enableShaderProcs(); // re-fetch ptrs in new context.
  }

  SDL_Window *window() { return window_; }
  int width() const { return width_; }
  int height() const { return height_; }

private:
  int display_;
  int width_, height_;           // current dimensions, window or fullscreen.
  int last_x_, last_y_;          // last known position of window
  int last_width_, last_height_; // last known dimension of window.
  bool fullscreen_;
  int ndisplays_;
  SDL_Window *window_;
  SDL_GLContext glcontext_;
  SDL_DisplayMode mode_[kMAXDISPLAYS];
  SDL_Rect rect_[kMAXDISPLAYS];
} window;

// Optional #defines for glsl compilation from .cfg file.
string defines;

StereoMode stereoMode = ST_NONE;

Uniforms uniforms;

string lifeform; // Conway's Game of Life creature, if any.

// Render globals grouped.
#define NFBO 2
#define NBLUR 2

struct RenderContext {
  GLuint mainFbo[NFBO];
  GLuint mainTex[NFBO];
  GLuint mainDepth[NFBO];

  GLuint fxaaFbo = -1;
  GLuint fxaaTex = -1;

  GLuint scratchFbo = -1;
  GLuint scratchTex = -1;

  GLuint blurFbo[NBLUR];
  GLuint blurTex[NBLUR];

  ShaderManager shaderManager;

  GLuint de_fbo;
  GLuint de_texture;

  GLuint background_texture;

  // DLP-Link or interlaced L/R eye polarity
  int polarity = 1;
} render;

// Try release everything that might have been allocated within
// current glContext. Leaks detected with AMD CodeXL.
void clearGlContext() {
  // release shaders.
  ::render.shaderManager.fractal.clear();
  ::render.shaderManager.effects.clear();
  ::render.shaderManager.dof.clear();
  ::render.shaderManager.de_shader.clear();
  ::render.shaderManager.fxaa.clear();

  // delete fbos and textures.
  glDeleteTextures(1, &render.de_texture);
  glDeleteFramebuffers(1, &render.de_fbo);

  glDeleteTextures(1, &render.background_texture); // free existing

  glDeleteFramebuffers(ARRAYSIZE(render.mainFbo), render.mainFbo);
  glDeleteTextures(ARRAYSIZE(render.mainDepth), render.mainDepth);
  glDeleteTextures(ARRAYSIZE(render.mainTex), render.mainTex);

  glDeleteFramebuffers(ARRAYSIZE(render.blurFbo), render.blurFbo);
  glDeleteTextures(ARRAYSIZE(render.blurTex), render.blurTex);

  glDeleteFramebuffers(1, &render.scratchFbo);
  glDeleteTextures(1, &render.scratchTex);

  glDeleteFramebuffers(1, &render.fxaaFbo);
  glDeleteTextures(1, &render.fxaaTex);
}

struct InputContext {
  bool grabbed = false;
  SDL_Joystick *stick = NULL;
  SDL_Cursor *arrow = NULL;
  SDL_Cursor *hand = NULL;
  SDL_Cursor *crosshair = NULL;
} input;

////////////////////////////////////////////////////////////////
// FPS tracking.
struct FPSCounter {
  int framesToAverage;
  std::vector<Uint32> frameDurations;
  int frameDurationsIndex = 0;
  Uint32 lastFrameTime;

  double now() { return (double)SDL_GetTicks() / 1000.0; }

  // Initialize the FPS structure.
  void init(int framesToAverage_) {
    assert(framesToAverage_ > 1);
    framesToAverage = framesToAverage_;
    frameDurations.resize(framesToAverage_, 0);
    lastFrameTime = SDL_GetTicks();
  }

  // Update the FPS structure after drawing a frame.
  void update(void) {
    Uint32 time = SDL_GetTicks();
    frameDurations[frameDurationsIndex++ % framesToAverage] =
        time - lastFrameTime;
    lastFrameTime = time;
  }

  // Return the duration of the last frame.
  Uint32 getLastFrameDuration(void) {
    return frameDurations[(frameDurationsIndex + framesToAverage - 1) %
                          framesToAverage];
  }

  // Return the average FPS over the last X frames.
  float get() {
    if (frameDurationsIndex < framesToAverage)
      return 0; // not enough data
    Uint32 sum = 0;
    for (int i = 0; i < framesToAverage; i++)
      sum += frameDurations[i];
    float fps = framesToAverage * 1000.f / sum;

    static Uint32 lastfps = 0;
    if (lastFrameTime - lastfps > 5000) {
      // Once per 5 seconds.
#if defined(_WIN32)
      // SetOculusPrediction(0.9 / fps); // A bit lower than latency.
#endif
      printf("fps %f\n", fps);
      lastfps = lastFrameTime;
    }

    return fps;
  }
} fpsCounter;

////////////////////////////////////////////////////////////////
// Current logical state of the program.

Camera camera, // Active camera view.
    config;    // Global configuration set.

vector<KeyFrame> keyframes; // Keyframes

void suggestDeltaTime(KeyFrame &camera, const vector<KeyFrame> &keyframes) {
  if (keyframes.empty()) {
    camera.delta_time = 0;
  } else {
    double dist = camera.distanceTo(keyframes[keyframes.size() - 1]);
    double steps = dist / camera.speed;
    camera.delta_time = steps / config.fps;
  }
}

#define NSUBFRAMES 100 // # splined subframes between keyframes.
// TODO: make relative to delta_time, thus more like max fps.

void CatmullRom(const vector<KeyFrame> &keyframes, vector<KeyFrame> *output,
                bool loop = false, int nsubframes = NSUBFRAMES) {
  output->clear();
  if (keyframes.size() < 2)
    return; // Need at least two points.

  vector<KeyFrame> controlpoints(keyframes);

  size_t n = controlpoints.size();

  // Compute / check quats.
  // Pick smaller angle between two quats.
  mat2quat(controlpoints[0].v, controlpoints[0].q);
  quat2x(controlpoints[0].q, controlpoints[0].x);
  for (size_t i = 1; i < n; ++i) {
    mat2quat(controlpoints[i].v, controlpoints[i].q);
    double dot = controlpoints[i - 1].q[0] * controlpoints[i].q[0] +
                 controlpoints[i - 1].q[1] * controlpoints[i].q[1] +
                 controlpoints[i - 1].q[2] * controlpoints[i].q[2] +
                 controlpoints[i - 1].q[3] * controlpoints[i].q[3];
    if (dot < 0) {
      // Angle between quats > 180; make current quat go smaller 360 - angle.
      // Note: this fails to cover the two ends of a loop.
      controlpoints[i].q[0] *= -1;
      controlpoints[i].q[1] *= -1;
      controlpoints[i].q[2] *= -1;
      controlpoints[i].q[3] *= -1;
    }
    // Compute splinable R4 mapping of quat for splining below.
    quat2x(controlpoints[i].q, controlpoints[i].x);
  }

  if (loop) {
    // Replicate first two at end to function as p2, p3 for splining.
    controlpoints.push_back(controlpoints[0]);
    if (controlpoints[controlpoints.size() - 1].delta_time == 0) {
      // Likely first point does not have delta_time specified.
      // Try guess at one based on speed at last point and distance there to.
      suggestDeltaTime(controlpoints[controlpoints.size() - 1], keyframes);
    }
    controlpoints.push_back(controlpoints[1]);
    // Last specified keyframe is p0 for spline of first keyframe.
    controlpoints.push_back(controlpoints[n - 1]);
  } else {
    // Replicate last one twice more.
    controlpoints.push_back(controlpoints[n - 1]);
    controlpoints.push_back(controlpoints[n - 1]);
    // Last one is p0 for spline of first keyframe.
    controlpoints.push_back(controlpoints[0]);
  }

  n = controlpoints.size();

  // Compute each frame's target time based on sum of delta_time up to it.
  // Note we don't spline delta_time but we do spline time.
  double time = 0;
  controlpoints[0].time = 0; // time starts at 0.
  for (size_t i = 1; i < n - 1; ++i) {
    time += controlpoints[i].delta_time;
    controlpoints[i].time = time;
  }
  // Last one's time is p0 for spline of first one; set to 0 as well.
  controlpoints[n - 1].time = 0;

  // Now spline all into intermediate frames.
  for (size_t i = 0; i < n - 3; ++i) {
    const KeyFrame *p0 = &controlpoints[i > 0 ? i - 1 : n - 1];
    const KeyFrame *p1 = &controlpoints[i];
    const KeyFrame *p2 = &controlpoints[i + 1];
    const KeyFrame *p3 = &controlpoints[i + 2];
    for (int f = 0; f < nsubframes; ++f) {
      KeyFrame tmp = config; // Start with default values.
      tmp.setKey(f == 0);
      const double t = ((double)f) / nsubframes;

// The CatmullRom spline function; 0 <= t <= 1
// Suffers from overshoot for non-evenly spaced control points.
// TODO: look into Bessel-Overhauser mitigation.
#define SPLINE(X, p0, p1, p2, p3)                                              \
  ((X) = (double)(.5 * (2 * (p1) +                                             \
                        t * ((-(p0) + (p2)) +                                  \
                             t * ((2 * (p0)-5 * (p1) + 4 * (p2) - (p3)) +      \
                                  t * (-(p0) + 3 * (p1)-3 * (p2) + (p3)))))))

      // Spline over splinable representation of quat.
      for (size_t j = 0; j < 4; ++j) {
        SPLINE(tmp.x[j], p0->x[j], p1->x[j], p2->x[j], p3->x[j]);
      }
      x2quat(tmp.x, tmp.q); // convert back to quat
      qnormalize(tmp.q);
      quat2mat(tmp.q, tmp.v); // convert quat to the splined rotation matrix

      // Spline position into tmp.v[12..14]
      for (size_t j = 12; j < 15; ++j) {
        // To control numerical precision, re-base to (p2-p1)/2.
        double a = p0->v[j], b = p1->v[j], c = p2->v[j], d = p3->v[j];
        double base = .5 * (c - b);
        a -= base;
        b -= base;
        c -= base;
        d -= base;
        SPLINE(tmp.v[j], a, b, c, d);
        tmp.v[j] += base;
      }

      // Spline par[] array. Some of those could also be rotations,
      // which will not spline nicely at all times..
      // TODO: have couple of uniform quats for shader use and spline those.
      for (size_t j = 0; j < ARRAYSIZE(tmp.par); ++j) {
        SPLINE(tmp.par[j][0], p0->par[j][0], p1->par[j][0], p2->par[j][0],
               p3->par[j][0]);
        SPLINE(tmp.par[j][1], p0->par[j][1], p1->par[j][1], p2->par[j][1],
               p3->par[j][1]);
        SPLINE(tmp.par[j][2], p0->par[j][2], p1->par[j][2], p2->par[j][2],
               p3->par[j][2]);
      }

      // Spline generic float uniform array.
      for (int j = 0; j < tmp.n_funis; ++j) {
        SPLINE(tmp.funis[j], p0->funis[j], p1->funis[j], p2->funis[j],
               p3->funis[j]);
      }

      // Spline common params, if marked as such.
#define PROCESS(a, b, c, doSpline)                                             \
  if (doSpline) {                                                              \
    SPLINE(tmp.b, p0->b, p1->b, p2->b, p3->b);                                 \
  }
      PROCESS_COMMON_PARAMS
#undef PROCESS

#undef SPLINE

      tmp.orthogonalize();
      output->push_back(tmp);
    }
  }
}

////////////////////////////////////////////////////////////////
// Controllers.

typedef enum Controller {
  // user parameters: 0..9
  // other parameters:
  CTL_FOV = ARRAYSIZE(camera.par),
  CTL_RAY,
  CTL_ITER,
  CTL_AO,
  CTL_GLOW,
  CTL_TIME,
  CTL_CAM,
  CTL_3D,
  CTL_LAST = CTL_3D,
} Controller;

// Controller modifiers.
// They get a pointer to the modified value and a signed count of consecutive
// changes.
void m_mul(float *x, int d) { *x *= pow(10, GLSL::sign(d) / 20.); }
void m_mulSlow(float *x, int d) { *x *= pow(10, GLSL::sign(d) / 40.); }
void m_mulSlow(double *x, int d) { *x *= pow(10, GLSL::sign(d) / 40.); }
void m_tan(float *x, int d) {
  *x = atan(tan(*x * PI / 180 / 2) * pow(0.1, GLSL::sign(d) / 40.)) / PI * 180 *
       2;
}
void m_progressiveInc(int *x, int d) {
  *x += GLSL::sign(d) * ((abs(d) + 4) / 4);
}
void m_progressiveAdd(float *x, int d) {
  *x += 0.001 * (GLSL::sign(d) * ((abs(d) + 4) / 4));
}
void m_progressiveAdd(double *x, int d) {
  *x += 0.001 * (GLSL::sign(d) * ((abs(d) + 4) / 4));
}
void m_singlePress(int *x, int d) {
  if (d == 1 || d == -1)
    *x += d;
}

void m_rotateX(int d) {
  camera.rotate(GLSL::sign(d) * camera.keyb_rot_speed, 0, 1, 0);
}
void m_rotateY(int d) {
  camera.rotate(-GLSL::sign(d) * camera.keyb_rot_speed, 1, 0, 0);
}

void m_rotateX2(float d) { camera.rotate(d, 0, 1, 0); }
void m_rotateY2(float d) { camera.rotate(d, 1, 0, 0); }
void m_rotateZ2(float d) { camera.rotate(d, 0, 0, 1); }

// Print controller values into a string.
char *printController(char *s, Controller c) {
  assert(c <= CTL_LAST);
  switch (c) {
  default: {
    char x[20], y[20];
    sprintf(x, "par%dx", c);
    sprintf(y, "par%dy", c);
    sprintf(s, "%s %.3f %s %.3f", y, camera.par[c][1], x, camera.par[c][0]);
  } break;
  case CTL_FOV:
    sprintf(s, "Fov %.3g %.3g", camera.fov_x, camera.fov_y);
    break;
  case CTL_RAY:
    sprintf(s, "Ray %.2e steps %d", camera.min_dist, camera.max_steps);
    break;
  case CTL_ITER:
    sprintf(s, "It %d|%d", camera.iters, camera.color_iters);
    break;
  case CTL_AO:
    sprintf(s, "aO %.2e aOeps %.2e", camera.ao_strength, camera.ao_eps);
    break;
  case CTL_GLOW:
    sprintf(s, "Glow %.3f Dist %.2e", camera.glow_strength,
            camera.dist_to_color);
    break;
  case CTL_CAM: {
    sprintf(s, "Look [%4d %4d %4d]", (int)(camera.ahead()[0] * 100),
            (int)(camera.ahead()[1] * 100), (int)(camera.ahead()[2] * 100));
  } break;
  case CTL_TIME: {
    sprintf(s, "Speed %.8e DeltaT %.3f", camera.speed, camera.delta_time);
  } break;
  case CTL_3D: {
    sprintf(s, "Sep %.8e Foc %.3f", camera.speed, camera.focus);
  } break;
  }
  return s;
}

// Update controller.y by the signed count of consecutive changes.
void updateControllerY(Controller c, int d, bool alt) {
  assert(c <= CTL_LAST);
  switch (c) {
  default:
    m_progressiveAdd(&camera.par[c][1], d);
    break;
  case CTL_FOV:
    m_tan(&camera.fov_x, d);
    m_tan(&camera.fov_y, d);
    break;
  case CTL_RAY:
    m_mul(&camera.min_dist, d);
    break;
  case CTL_ITER:
    m_singlePress(&camera.iters, d);
    if ((camera.iters & 1) == (d > 0))
      m_singlePress(&camera.color_iters, d);
    break;
  case CTL_AO:
    m_mulSlow(&camera.ao_strength, d);
    break;
  case CTL_GLOW:
    m_progressiveAdd(&camera.glow_strength, d);
    break;
  case CTL_CAM:
    m_rotateY(d);
    break;
  case CTL_TIME:
    m_progressiveAdd(&camera.delta_time, d);
    break;
  case CTL_3D:
    m_progressiveAdd(&camera.focus, d * 100);
    break;
  }
  // Enforce sane bounds.
  if (camera.delta_time < 0)
    camera.delta_time = 0;
}

// Update controller.x by the signed count of consecutive changes.
void updateControllerX(Controller c, int d, bool alt) {
  assert(c <= CTL_LAST);
  switch (c) {
  default:
    m_progressiveAdd(&camera.par[c][alt ? 2 : 0], d);
    break;
  case CTL_FOV:
    m_tan(&camera.fov_x, d);
    break;
  case CTL_RAY:
    m_progressiveInc(&camera.max_steps, d);
    break;
  case CTL_ITER:
    m_singlePress(&camera.color_iters, d);
    break;
  case CTL_AO:
    m_mul(&camera.ao_eps, d);
    break;
  case CTL_GLOW:
    m_mul(&camera.dist_to_color, d);
    break;
  case CTL_CAM:
    m_rotateX(d);
    break;
  case CTL_TIME:
    m_mulSlow(&camera.speed, d);
    break;
  case CTL_3D:
    m_mulSlow(&camera.speed, d);
    break;
  }
}

// Change the active controller with a keypress.
void changeController(SDL_Keycode key, Controller *c) {
  switch (key) {
  case SDLK_0:
  case SDLK_KP_0:
    *c = (Controller)0;
    break;
  case SDLK_1:
  case SDLK_KP_1:
    *c = (Controller)1;
    break;
  case SDLK_2:
  case SDLK_KP_2:
    *c = (Controller)2;
    break;
  case SDLK_3:
  case SDLK_KP_3:
    *c = (Controller)3;
    break;
  case SDLK_4:
  case SDLK_KP_4:
    *c = (Controller)4;
    break;
  case SDLK_5:
  case SDLK_KP_5:
    *c = (Controller)5;
    break;
  case SDLK_6:
  case SDLK_KP_6:
    *c = (Controller)6;
    break;
  case SDLK_7:
  case SDLK_KP_7:
    *c = (Controller)7;
    break;
  case SDLK_8:
  case SDLK_KP_8:
    *c = (Controller)8;
    break;
  case SDLK_9:
  case SDLK_KP_9:
    *c = (Controller)9;
    break;
  case SDLK_f:
    *c = CTL_FOV;
    break;
  case SDLK_g:
    *c = CTL_GLOW;
    break;
  case SDLK_l:
    *c = CTL_CAM;
    break;
  case SDLK_i:
    *c = CTL_ITER;
    break;
  case SDLK_o:
    *c = CTL_AO;
    break;
  case SDLK_r:
    *c = CTL_RAY;
    break;
  case SDLK_t:
    *c = CTL_TIME;
    break;
  case SDLK_v:
    *c = CTL_3D;
    break;
  default:
    break; // no change
  }
}

////////////////////////////////////////////////////////////////
// Graphics.

void saveScreenshot(char const *tgaFile) {
  TGA tga;
  tga.fromFramebuffer(config.width, config.height);
  string filename(WorkingDir + tgaFile);
  if (tga.writeFile(filename.c_str()))
    printf(__FUNCTION__ " : wrote %s\n", filename.c_str());
  else
    printf(__FUNCTION__ " : failed to write %s\n", filename.c_str());
}

TGA background;

void LoadBackground() {
  string filename(WorkingDir + "background.tga");
  background.readFile(filename.c_str());
  if (background.data()) {
    printf(__FUNCTION__ " : loaded background image from '%s'\n",
           filename.c_str());
  }
}

// Return BGR value of pixel x,y.
unsigned int getBGRpixel(int x, int y) {
  unsigned char img[3];
  int height = config.height;
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_FRONT);
  glReadPixels(x, height - 1 - y, 1, 1, GL_BGR, GL_UNSIGNED_BYTE, img);
  unsigned int val = img[0] * 256 * 256 + img[1] * 256 + img[2];
  return val;
}

GLSL::vec3 getPixelColor(int x, int y) {
  float img[3];
  int height = config.height;
  glReadBuffer(GL_COLOR_ATTACHMENT0);
  glReadPixels(x, height - 1 - y, 1, 1, GL_RGB, GL_FLOAT, img);
  return GLSL::vec3(img);
}

bool setupDirectories(const char *configFile) {
  BaseFile.clear();
  BaseDir.clear();
  WorkingDir.clear();

  char dirName[MAX_PATH];
#if defined(_WIN32)
  GetModuleFileName(NULL, dirName, MAX_PATH);
  BaseDir.assign(dirName);
  size_t slash = BaseDir.rfind('\\');
  if (slash != string::npos)
    BaseDir.erase(slash + 1);

  char *fileName = NULL;
  DWORD result = GetFullPathName(configFile, MAX_PATH, dirName, &fileName);
  if (result) {
    if (fileName) {
      BaseFile.assign(fileName);
      *fileName = '\0';
    }
    WorkingDir.assign(dirName);
  }
#else // !_WIN32
  strncpy(dirName, configFile, sizeof dirName);
  dirName[sizeof dirName - 1] = 0;
  int i = strlen(dirName);
  while (i > 0 && strchr("/", dirName[i - 1]) == NULL)
    --i;
  BaseFile.assign(&dirName[i]);
  dirName[i] = 0;
  BaseDir.assign("./"); // Assume relative to cwd.
  WorkingDir.assign(dirName);
#endif

  if (BaseFile.empty())
    BaseFile.assign(DEFAULT_CONFIG_FILE);

  cout << __FUNCTION__ << " : " << BaseDir << ", " << WorkingDir << ", "
       << BaseFile << endl;

  return true;
}

// Initializes the video mode, OpenGL state, shaders,
// camera and shader parameters.
// Grabs mouse input.
//
// Exits the program if an error occurs.
bool initGraphics(bool fullscreenToggle, int w, int h, bool hideMouse) {
  // Set attributes for the OpenGL window.
  SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, config.depth_size);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

  if (stereoMode == ST_QUADBUFFER) {
    SDL_GL_SetAttribute(SDL_GL_STEREO, 1);
  }

  int ww = w, hh = h;
  if (w == h && (h == 2048 || h == 4096)) {
    ww = 1024; // HACK!
    hh = 1024;
  }
  if (w == 2 * h && (w == 4096 || w == 2048)) {
    ww = 2048; // HACK!
    hh = 1024;
  }

  if (fullscreenToggle) {
    if (!window.toggleFullscreen()) {
      return false;
    }
  } else {
    if (!window.resize(ww, hh)) {
      // Spurious SDL resize? Ignore.
      return false;
    }
  }

  GLint maxRenderBufferSize;
  glGetIntegerv(GL_MAX_RENDERBUFFER_SIZE, &maxRenderBufferSize);
  printf(__FUNCTION__ " : max render buffer size %d\n", maxRenderBufferSize);
  GLint dims[2];
  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, &dims[0]);
  printf(__FUNCTION__ " : max viewport dims %dx%d\n", dims[0], dims[1]);

  printf(__FUNCTION__ " : GetSwap %d\n", SDL_GL_GetSwapInterval());

  if (stereoMode == ST_QUADBUFFER) {
    int ga = 0;
    SDL_GL_GetAttribute(SDL_GL_STEREO, &ga);
    if (ga == 0)
      die("No stereo rendering available: %s\n", SDL_GetError());
  }

  // Take the window size the system provided. Might be higher than requested.
  if (ww == w && hh == h) {
    config.width = window.width();
    config.height = window.height();
  }

  printf(__FUNCTION__ " : %dx%d\n", config.width, config.height);

  SDL_GL_GetAttribute(SDL_GL_DEPTH_SIZE, &config.depth_size);
  printf(__FUNCTION__ " : depth size %u\n", config.depth_size);

  if (input.arrow == NULL)
    input.arrow = SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_ARROW);
  if (input.hand == NULL)
    input.hand = SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_HAND);
  if (input.crosshair == NULL)
    input.crosshair = SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_CROSSHAIR);

  input.grabbed = true;

  if (hideMouse) {
    XXX_SDL_SetCursor(input.arrow);
    SDL_SetRelativeMouseMode(SDL_FALSE); // This toggle is needed!
    SDL_SetRelativeMouseMode(SDL_TRUE);
  } else {
    XXX_SDL_SetCursor(input.crosshair);
  }

  printf(__FUNCTION__ " : SetSwap %d\n", SDL_GL_SetSwapInterval(0));

  // Enable shader functions and compile shaders.
  // Needs to be done after setting the video mode.
  enableShaderProcs() || die("This program needs support for GLSL shaders.\n");

  ::render.shaderManager.loadFractal(defines) ||
      die("Error in GLSL ::render.shaderManager.fractal shader "
          "compilation:\n%s\n",
          ::render.shaderManager.fractal.log().c_str());

  ::render.shaderManager.loadHelpers(defines, stereoMode) ||
      die("Error in GLSL ::render.shaderManager.effects shader "
          "compilation:\n%s\n",
          ::render.shaderManager.effects.log().c_str());

  glEnable(GL_TEXTURE_2D);

  if (!config.disable_de) {
#if defined(GL_RGBA32F) // We need this to be actually capable of GL_FLOAT
    // Try compile same shader to get a minimal DE computation version.
    ::render.shaderManager.loadDEShader(defines,
                                        "#define ST_COMPUTE_DE_ONLY\n");

    if (!::render.shaderManager.de_shader.ok()) {
      printf(__FUNCTION__ " : ::render.shaderManager.de_shader failed to "
                          "compile: no GPU de.\n");
      printf(__FUNCTION__ " :\n%s\n",
             ::render.shaderManager.de_shader.log().c_str());
      while (glGetError() != GL_NO_ERROR)
        ;
    } else {
      // Got a shader that can compute DE.
      // Set-up a float32 fbo for it to write to.
      glGenTextures(1, &render.de_texture);

      glGenFramebuffers(1, &render.de_fbo);

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, render.de_texture);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, config.width, config.height, 0,
                   GL_RGB, GL_FLOAT, NULL);
      glBindTexture(GL_TEXTURE_2D, 0);

      glBindFramebuffer(GL_FRAMEBUFFER, render.de_fbo);

      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                             GL_TEXTURE_2D, render.de_texture, 0);

      CHECK_ERROR;
      CHECK_FRAMEBUFFER;

      glClear(GL_COLOR_BUFFER_BIT);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
#endif // defined(RGBA32F)
  }

// Play hokey w/ substandard osx ogl support.
#if defined(GL_RGBA16F)
  // try half-float first. Faster and enough range?
#define xGL_RGBA32F GL_RGBA16F
#define xGL_FLOAT GL_HALF_FLOAT
#elif defined(GL_RGBA32F)
  // try float next.
#define xGL_RGBA32F GL_RGBA32F
#define xGL_FLOAT GL_FLOAT
#else
  // no HDR fbo for you..
#define xGL_RGBA32F GL_RGBA16
#define xGL_FLOAT GL_UNSIGNED_SHORT
#endif

  if (background.data() != NULL) {

    // Load background image into texture.

    glGenTextures(1, &render.background_texture);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, render.background_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    GLint magFilter = config.backbuffer ? GL_NEAREST : GL_LINEAR;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8, // assume tga is srgb
                 background.width(), background.height(), 0, GL_BGR,
                 GL_UNSIGNED_BYTE, background.data());
    glGenerateMipmap(GL_TEXTURE_2D);
    printf(__FUNCTION__ " : background texture at %d\n",
           render.background_texture);
    glBindTexture(GL_TEXTURE_2D, 0);
    CHECK_ERROR;
  }

  if (config.enable_dof || config.backbuffer) {
    // Set up fbos for multipass rendering.

    // Create depth buffer(s)
    glGenTextures(ARRAYSIZE(render.mainDepth), render.mainDepth);

    // Create textures to render to
    glGenTextures(ARRAYSIZE(render.mainTex), render.mainTex);

    // Create framebuffers
    glGenFramebuffers(ARRAYSIZE(render.mainFbo), render.mainFbo);

    GLint clamp = config.backbuffer ? GL_REPEAT : GL_CLAMP_TO_EDGE;
    GLint minfilter = config.backbuffer ? GL_NEAREST : GL_LINEAR;

    for (size_t i = 0; i < ARRAYSIZE(render.mainFbo); ++i) {
      // Initialize depth texture
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, render.mainDepth[i]);

      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, config.width,
                   config.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

      CHECK_ERROR;

      // Initialize color texture
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, render.mainTex[i]);

      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, clamp);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, clamp);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minfilter);

      // 32F gives best fidelity for the pathtracer and radiance shaders that
      // use accumulation.
#if defined(GL_RGBA32F)
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, config.width, config.height, 0,
                   GL_BGRA, GL_FLOAT, NULL);
#else
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, config.width, config.height, 0,
                   GL_BGRA, GL_FLOAT, NULL);
#endif

      CHECK_ERROR;

      glBindTexture(GL_TEXTURE_2D, 0);

      // Initialize render buffer
      glBindFramebuffer(GL_FRAMEBUFFER, render.mainFbo[i]);

      // Attach colorbuffer
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                             GL_TEXTURE_2D, render.mainTex[i], 0);

      CHECK_ERROR;

      // Attach depthbuffer
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                             render.mainDepth[i], 0);

      CHECK_ERROR;
      CHECK_FRAMEBUFFER;

      // Clear it all.
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    if (::render.shaderManager.dof.ok()) {
      // Create and initialize blur fbos.

      glGenTextures(ARRAYSIZE(render.blurTex), render.blurTex);

      for (size_t i = 0; i < ARRAYSIZE(render.blurTex); ++i) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, render.blurTex[i]);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

        glTexImage2D(GL_TEXTURE_2D, 0, xGL_RGBA32F, config.width, config.height,
                     0, GL_BGRA, xGL_FLOAT, NULL);

        glBindTexture(GL_TEXTURE_2D, 0);

        CHECK_ERROR;
      }

      // Create framebuffers to render to blur array texture.
      glGenFramebuffers(ARRAYSIZE(render.blurFbo), render.blurFbo);

      for (size_t i = 0; i < ARRAYSIZE(render.blurFbo); ++i) {
        glBindFramebuffer(GL_FRAMEBUFFER, render.blurFbo[i]);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D, render.blurTex[i], 0);

        CHECK_FRAMEBUFFER;
        CHECK_ERROR;

        glClear(GL_COLOR_BUFFER_BIT);

        // Back to normal framebuffer.
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
      }
    }

    if (::render.shaderManager.dof.ok() || ::render.shaderManager.fxaa.ok()) {

      // Create and initialize scratch fbo.

      glGenTextures(1, &render.scratchTex);
      glGenFramebuffers(1, &render.scratchFbo);

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, render.scratchTex);

      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

      glTexImage2D(GL_TEXTURE_2D, 0, xGL_RGBA32F, config.width, config.height,
                   0, GL_BGRA, xGL_FLOAT, NULL);

      glBindFramebuffer(GL_FRAMEBUFFER, render.scratchFbo);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                             GL_TEXTURE_2D, render.scratchTex, 0);

      CHECK_FRAMEBUFFER;
      CHECK_ERROR;

      glClear(GL_COLOR_BUFFER_BIT);

      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glBindTexture(GL_TEXTURE_2D, 0);
    } // ::render.shaderManager.dof.ok() || ::render.shaderManager.fxaa.ok()

    if (::render.shaderManager.fxaa.ok()) {

      // Create and initialize ::render.shaderManager.fxaa fbo

      glGenTextures(1, &render.fxaaTex);
      glGenFramebuffers(1, &render.fxaaFbo);

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, render.fxaaTex);

      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

      glTexImage2D(GL_TEXTURE_2D, 0, xGL_RGBA32F, config.width, config.height,
                   0, GL_BGRA, xGL_FLOAT, NULL);

      glBindFramebuffer(GL_FRAMEBUFFER, render.fxaaFbo);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                             GL_TEXTURE_2D, render.fxaaTex, 0);

      CHECK_FRAMEBUFFER;
      CHECK_ERROR;

      glClear(GL_COLOR_BUFFER_BIT);

      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glBindTexture(GL_TEXTURE_2D, 0);
    } // ::render.shaderManager.fxaa.ok()
  }

  // Fill backbuffer w/ starting lifeform, if we have one.

  if (config.backbuffer && !lifeform.empty()) {
    glBindFramebuffer(GL_FRAMEBUFFER, render.mainFbo[1 /* backbuffer */]);
    // Ortho projection, entire screen in regular pixel coordinates.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, config.width, config.height, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glPointSize(1);
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glColor4f(1, 1, 1, 1);
    glBegin(GL_POINTS);

    {
      istringstream in(lifeform);
      string line;
      bool isRLE = false, done = !in.good();
      int x = 0, y = 0, accu = 0, xo = 0, yo = 0;
      while (!done) {
        if (!isRLE) {
          if (!getline(in, line)) {
            done = true;
            continue;
          }
          if (line.empty() || line[0] == '#')
            continue;
          if (line[0] == 'x') {
            // Try parse RLE format lifeform dimensions and center it.
            sscanf(line.c_str(), "x = %d, y = %d,", &xo, &yo);
            xo = (config.width - xo) / 2;
            yo = (config.height - yo) / 2;
            if (xo < 0)
              xo = 0;
            if (yo < 0)
              yo = 0;
            isRLE = true;
            continue;
          }
          // Life 0.6 format.
          istringstream is(line);
          is >> x;
          is >> y;
          glVertex2f(.5 + config.width / 2 + x, .5 + config.height / 2 + y);
        } else {
          // RLE format.
          char c = in.get();
          if (!in.good()) {
            done = true;
            continue;
          }
          switch (c) {
          case '0':
          case '1':
          case '2':
          case '3':
          case '4':
          case '5':
          case '6':
          case '7':
          case '8':
          case '9': {
            accu *= 10;
            accu += c - '0';
          }
            continue;
          case '$': {
            if (!accu)
              accu = 1;
            x = 0;
            y += accu;
            accu = 0;
          }
            continue;
          case 'b':
          case 'B': {
            if (!accu)
              accu = 1;
            x += accu;
            accu = 0;
          }
            continue;
          case 'o':
          case 'A': {
            if (!accu)
              accu = 1;
            for (int xx = x; xx < x + accu; ++xx) {
              glVertex2f(.5 + xx + xo, .5 + y + yo);
            }
            x += accu;
            accu = 0;
          }
            continue;
          case '!':
            done = true;
            continue;
          }
        }
      }
    }

    glEnd();
    glFinish();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  return true;
}

TwBar *bar = NULL;

// Find '\n#define foo par[x].z  // {twbar params}' in
// ::render.shaderManager.glsl_source.
void initTwParDefines() {
  size_t start = 0;
  while ((start = ::render.shaderManager.glsl_source.find(
              "\n#define ", start + 1)) != string::npos) {
    size_t eol = ::render.shaderManager.glsl_source.find("\n", start + 1);
    if (eol == string::npos)
      continue;
    string line(::render.shaderManager.glsl_source, start + 1,
                eol - (start + 1));

    size_t parStart = line.find(" par[");
    if (parStart == string::npos || parStart < 8)
      continue;
    string varName(line, 8, parStart - 8);

    size_t attr_start = line.find("{");
    size_t attr_end = line.find("}");
    string attr;
    if (attr_start != string::npos && attr_end != string::npos &&
        attr_end > attr_start)
      attr.assign(line, attr_start + 1, attr_end - (attr_start + 1));

    int index;
    char _xyz = 'x';
    if (sscanf(line.c_str() + parStart + 5, "%d].%c", &index, &_xyz) < 1)
      continue;
    if (index < 0 || index > (int)ARRAYSIZE(camera.par))
      continue;
    if (_xyz < 'x' || _xyz > 'z')
      continue;

    printf("parameter %s par[%d].%c {%s}\n", varName.c_str(), index, _xyz,
           attr.c_str());

    float *address = &camera.par[index][_xyz - 'x'];

    if (varName.find("Color") != string::npos) {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_COLOR3F, address, "");
    } else if (varName.find("Vector") != string::npos) {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_DIR3F, address, "");
    } else if (varName.find("Quat") != string::npos) {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_QUAT4F, address, "");
    } else {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_FLOAT, address, attr.c_str());
    }
  }
}

// Find '\nuniform <type> name;' in ::render.shaderManager.glsl_source.
// Pick up type and twbar params, if found.
void initTwUniform(const string &name, void *addr) {
  size_t start = 0;
  while ((start = ::render.shaderManager.glsl_source.find(
              "\nuniform ", start + 1)) != string::npos) {
    size_t eol = ::render.shaderManager.glsl_source.find("\n", start + 1);
    string line(::render.shaderManager.glsl_source, start + 1,
                eol - (start + 1));

    size_t attr_start = line.find("{");
    size_t attr_end = line.find("}");
    string attr;
    if (attr_start != string::npos && attr_end != string::npos &&
        attr_end > attr_start)
      attr.assign(line, attr_start + 1, attr_end - (attr_start + 1));

    if (line.find("uniform float " + name + ";") == 0) {
      TwAddVarRW(bar, name.c_str(), TW_TYPE_FLOAT, (float *)addr, attr.c_str());
    } else if (line.find("uniform int " + name + ";") == 0) {
      TwAddVarRW(bar, name.c_str(), TW_TYPE_INT32, (int *)addr, attr.c_str());
    }
  }
}

void initTwBar(enum StereoMode stereoMode) {
  if (bar == NULL)
    TwInit(TW_OPENGL, NULL);

  TwWindowSize(window.width(), window.height());

  if (bar != NULL)
    return;

  bar = TwNewBar("boxplorer");

  if (stereoMode == ST_OCULUS) {
    // Position HUD center for left eye.
    // TODO: make float
    char pos[100];
    int x = window.width();
    int y = window.height();
    sprintf(pos, "boxplorer position='%d %d'", x / 6, y / 4);
    TwDefine(pos);
    TwDefine("GLOBAL fontsize=3");
  }

  if (::render.shaderManager.fxaa.ok()) {
    TwAddVarRW(bar, "fxaa", TW_TYPE_BOOL32, &camera.fxaa, "group=post");
  }
  if (::render.shaderManager.dof.ok()) {
    TwAddVarRW(bar, "dof", TW_TYPE_BOOL32, &camera.enable_dof, "group=post");
    TwAddVarRW(bar, "aperture", TW_TYPE_FLOAT, &camera.aperture,
               "min=0.0 max=10.0 step=0.01 group=post");
  }
  if (::render.shaderManager.effects.ok() || ::render.shaderManager.fxaa.ok() ||
      ::render.shaderManager.dof.ok()) {
    // Global, thus on config, not camera.
    TwAddVarRW(bar, "exposure", TW_TYPE_FLOAT, &config.exposure,
               "min=0.0 max=5.0 step=0.01 group=post");
    TwAddVarRW(bar, "maxBright", TW_TYPE_FLOAT, &config.maxBright,
               "min=0.0 max=5.0 step=0.01 group=post");
    TwAddVarRW(bar, "gamma", TW_TYPE_FLOAT, &config.gamma,
               "min=0.0 max=5.0 step=0.01 group=post");
  }

  TwAddVarRW(bar, "focus", TW_TYPE_FLOAT, &camera.focus,
             "min=-20.0 max=30.0 step=0.1 group=3d");

  if (stereoMode == ST_OCULUS) {
    TwAddVarRW(bar, "ipd", TW_TYPE_FLOAT, &config.ipd,
               "min=-10.0 max=10.0 step=0.1 group=oculus");
  }

  uniforms.bindToUI(bar);

  initTwParDefines();

  // Tweak menus depending on state.
  TwDefine("boxplorer/post opened=false");
  TwDefine("boxplorer/3d opened=false");
  TwDefine("boxplorer/oculus opened=false");

  if (!config.enable_dof) {
    TwDefine("boxplorer/post visible=false");
  }
  if (stereoMode != ST_OCULUS) {
    TwDefine("boxplorer/oculus visible=false");
  }
}

void LoadKeyFrames(bool fixedFov) {
  char filename[256];
  for (int i = 0;; ++i) {
    sprintf(filename, "%s-%u.cfg", kKEYFRAME, i);
    // We load into global camera since that's where
    // the uniforms are bound to.
    if (!camera.loadConfig(filename))
      break;
    if (fixedFov) {
      camera.width = config.width;
      camera.height = config.height;
      camera.fov_x = config.fov_x;
      camera.fov_y = config.fov_y;
    }
    keyframes.push_back(camera);
  }
  printf(__FUNCTION__ " : loaded %lu keyframes\n",
         (unsigned long)keyframes.size());
}

void drawScreen(int leftRight = 3) {
  // Draw our texture covering entire screen, running the frame shader.
  if (leftRight & 1) {
    glBegin(GL_QUADS); // left
    glTexCoord2f(0, 1);
    glVertex2f(0, 0);
    glTexCoord2f(0, 0);
    glVertex2f(0, config.height);
    glTexCoord2f(.5, 0);
    glVertex2f(config.width / 2, config.height);
    glTexCoord2f(.5, 1);
    glVertex2f(config.width / 2, 0);
    glEnd();
  }
  if (leftRight & 2) {
    glBegin(GL_QUADS); // right
    glTexCoord2f(0.5, 1);
    glVertex2f(config.width / 2, 0);
    glTexCoord2f(0.5, 0);
    glVertex2f(config.width / 2, config.height);
    glTexCoord2f(1, 0);
    glVertex2f(config.width, config.height);
    glTexCoord2f(1, 1);
    glVertex2f(config.width, 0);
    glEnd();
  }
}

////////////////////////////////////////////////////////////////
// Setup, input handling and drawing.

int main(int argc, char **argv) {
  bool rendering = false;
  bool rendercubes = false;
  bool loop = false;
  bool useTime = false;
  bool configSpeed = false;
  bool fixedFov = false;
  bool fullscreen = false;
  char *lifeform_file = NULL;
  int enableDoF = 0;
  int disableDE = 0;
  int disableSpline = 0;
  bool xbox360 = true; // try find xbox360 controller
  int kJOYSTICK = 0;   // input.stick by index
  char *outputFilename = NULL;

  // Peel known options off the back..
  char *prev_arg = NULL;
  while (argc > 1) {
    if (!strcmp(argv[argc - 1], "--overunder")) {
      stereoMode = ST_OVERUNDER;
    } else if (!strcmp(argv[argc - 1], "--interlaced")) {
      stereoMode = ST_INTERLACED;
      defines.append("#define ST_INTERLACED\n");
    } else if (!strcmp(argv[argc - 1], "--xeyed")) {
      stereoMode = ST_XEYED;
    } else if (!strcmp(argv[argc - 1], "--sidebyside")) {
      stereoMode = ST_SIDEBYSIDE;
    } else if (!strcmp(argv[argc - 1], "--quadbuffer")) {
      stereoMode = ST_QUADBUFFER;
    } else if (!strcmp(argv[argc - 1], "--anaglyph")) {
      stereoMode = ST_ANAGLYPH;
      defines.append("#define ST_ANAGLYPH\n");
    } else if (!strcmp(argv[argc - 1], "--spherical")) {
      stereoMode = ST_SPHERICAL;
      defines.append("#define ST_SPHERICAL\n");
    } else if (!strcmp(argv[argc - 1], "--dome")) {
      stereoMode = ST_DOME;
      defines.append("#define ST_DOME\n");
    } else if (!strcmp(argv[argc - 1], "--oculus")) {
      stereoMode = ST_OCULUS;
      defines.append("#define ST_OCULUS\n");
    } else if (!strcmp(argv[argc - 1], "--render")) {
      rendering = true;
    } else if (!strcmp(argv[argc - 1], "--cubes")) {
      rendercubes = true;
    } else if (!strcmp(argv[argc - 1], "--time")) {
      useTime = true;
    } else if (!strcmp(argv[argc - 1], "--fullscreen")) {
      fullscreen = true;
    } else if (!strcmp(argv[argc - 1], "--speed")) {
      configSpeed = true;
    } else if (!strcmp(argv[argc - 1],
                       "--disable-::render.shaderManager.dof")) {
      enableDoF = -1;
    } else if (!strcmp(argv[argc - 1], "--enable-::render.shaderManager.dof")) {
      enableDoF = 1;
    } else if (!strcmp(argv[argc - 1], "--disable-de")) {
      disableDE = -1;
    } else if (!strcmp(argv[argc - 1], "--enable-de")) {
      disableDE = 1;
    } else if (!strcmp(argv[argc - 1], "--disable-spline")) {
      disableSpline = -1;
    } else if (!strcmp(argv[argc - 1], "--enable-spline")) {
      disableSpline = 1;
    } else if (!strcmp(argv[argc - 1], "--fixedfov")) {
      fixedFov = true;
    } else if (!strcmp(argv[argc - 1], "--loop")) {
      loop = true;
    } else if (!strcmp(argv[argc - 1], "--no360")) {
      xbox360 = false;
    } else if (!strncmp(argv[argc - 1], "--output=", 9)) {
      outputFilename = argv[argc - 1] + 9;
    } else if (!strncmp(argv[argc - 1], "--kf=", 5)) {
      kKEYFRAME = argv[argc - 1] + 5;
    } else if (!strncmp(argv[argc - 1], "--input.stick=", 11)) {
      kJOYSTICK = atoi(argv[argc - 1] + 11);
    } else if (!strncmp(argv[argc - 1], "--lifeform=", 11)) {
      lifeform_file = argv[argc - 1] + 11;
    } else if (!strcmp(argv[argc - 1], "--lifeform") && prev_arg != NULL) {
      lifeform_file = prev_arg;
    } else {
      prev_arg = argv[argc - 1];
      if (argc < 3)
        break;
    }
    --argc;
  }

#if defined(_WIN32)
  SetProcessDPIAware();
#endif

  if (stereoMode == ST_NONE)
    defines.append("#define ST_NONE\n");

  const char *configFile = (argc >= 2 ? argv[1] : DEFAULT_CONFIG);

  // Load configuration.
  if (setupDirectories(configFile) && config.loadConfig(BaseFile, &defines)) {
    // succuss
  } else {
    { die("Usage: boxplorer2 <configuration-file.cfg>\n"); }
  }

  if (lifeform_file) {
    // Try extract base name, since readFile loads from current data directory.
    // By allowed paths in the name cmdline completion can be used for the arg.
    const char *base = strrchr(lifeform_file, '/');
    if (!base)
      base = strrchr(lifeform_file, '\\');
    if (!base)
      base = lifeform_file;
    // Load definition into our global string.
    readFile(base, &lifeform);
  }

  if (fullscreen)
    config.fullscreen = fullscreen;

#if defined(_WIN32)
  if (stereoMode == ST_OCULUS) {
    if (!InitOculusSDK()) {
      die("InitOculusSDK() fail!");
    }

    hmd_settings_t hmd;
    if (!GetOculusDeviceInfo(&hmd)) {
      die("GetOculusDeviceInfo() fail!");
    }
    printf("Oculus %dx%d\n", hmd.h_resolution, hmd.v_resolution);
    printf("Oculus ipd %f\n", hmd.interpupillary_distance);
    printf("Oculus distortion (%f,%f,%f,%f)\n", hmd.distortion_k[0],
           hmd.distortion_k[1], hmd.distortion_k[2], hmd.distortion_k[3]);
    // TODO: actually do something w/ these values.

    SetOculusPrediction(.025); // Also gets adjusted later based on fps.
    ResetOculusOrientation();
  }
#endif

#if defined(HYDRA)
  if (sixenseInit() != SIXENSE_SUCCESS) {
    die("sixenseInit() fail!");
  }

  sixenseSetFilterEnabled(1);
  if (sixenseSetFilterParams(0.0, 0.0, 2000.0, 1.0) != SIXENSE_SUCCESS) {
    die("SetFilterParams() fail!");
  }

  if (sixenseSetActiveBase(0) != SIXENSE_SUCCESS) {
    die("sixenseSetActiveBase() fail!");
  }

  sixenseAllControllerData ssdata;

  if (sixenseIsControllerEnabled(0) != SIXENSE_SUCCESS) {
    die("controller(0) not enabled!");
  }

  double speed_base = 200.0;
  int lbuttons = SIXENSE_BUTTON_START; // so we calibrate on first loop.
  int rbuttons = 0;
  float neutral_x = 0;
  float neutral_y = 0;
  float neutral_z = 0;
#endif

  double speed_factor = 1.0;

  // Sanitize / override config parameters.
  if (loop)
    config.loop = true;
  if (enableDoF)
    config.enable_dof = (enableDoF == 1); // override
  if (stereoMode == ST_INTERLACED || stereoMode == ST_QUADBUFFER ||
      stereoMode == ST_ANAGLYPH) {
    config.enable_dof =
        0; // ::render.shaderManager.fxaa post does not work for these.
  }
  if (disableDE)
    config.disable_de = (disableDE == -1); // override
  if (disableSpline)
    config.no_spline = (disableSpline == -1); // override
  if (stereoMode == ST_OCULUS) {
    // Fix resolution for optimal performance.
    // config.width = 1280; config.height = 800;  // DK1
    // config.fov_x = 110; config.fov_y = 94.0;  // DK1
    config.width = 1920;
    config.height = 1080;
    config.fov_x = 100;
    config.fov_y = 100.0;
    fixedFov = true;
    // Enable multipass but not ::render.shaderManager.dof and
    // ::render.shaderManager.fxaa.
    config.backbuffer = 1;
    config.enable_fxaa = 0;
    config.enable_dof = 0;
  }
  if (stereoMode == ST_SPHERICAL) {
    config.width = 2048;
    config.height = config.width / 2;
  }
  if (stereoMode == ST_DOME) {
    config.height = config.width; // square
  }
  if (stereoMode == ST_INTERLACED || stereoMode == ST_OVERUNDER) {
    // Fix at 1080P
    config.width = 1920;
    config.height = 1080;
    config.fov_y = 30;
    config.fov_x = 0.0;
    fixedFov = true;
  }
  if (config.fps < 5)
    config.fps = 30;
  if (config.depth_size < 16)
    config.depth_size = 16;

  if (rendercubes) {
    // Render 6 cube faces:
    // make sure we have square view at 90 degrees.
    fixedFov = true;
    config.width = 2048;
    config.height = 2048;
    config.fov_x = 90.0;
    config.fov_y = 90.0;
    config.enable_fxaa = 1;
    config.enable_dof = 1;
    stereoMode = ST_NONE;
  }

  int saveWidth = config.width;
  int saveHeight = config.height;

  if (stereoMode == ST_XEYED) {
    // config.width *= 2;
    config.width = 3840;
    config.height = 1080;
  }

  config.sanitizeParameters();

  printf(__FUNCTION__ ": sanitized: size %dx%d\n", config.width, config.height);
  printf(__FUNCTION__ ": sanitized: view %fx%f\n", config.fov_x, config.fov_y);

  LoadBackground();

  // Initialize SDL and OpenGL graphics.
  SDL_Init(SDL_INIT_VIDEO) == 0 ||
      die("SDL initialization failed: %s\n", SDL_GetError());
  atexit(SDL_Quit);

  window.init();

  if (kJOYSTICK) {
    // open a input.stick by explicit index.
    SDL_InitSubSystem(SDL_INIT_JOYSTICK);
    input.stick = SDL_JoystickOpen(kJOYSTICK - 1);
    printf(__FUNCTION__ " : JoystickName '%s'\n",
           SDL_JoystickName(input.stick));
    printf(__FUNCTION__ " : JoystickNumAxes   : %i\n",
           SDL_JoystickNumAxes(input.stick));
    printf(__FUNCTION__ " : JoystickNumButtons: %i\n",
           SDL_JoystickNumButtons(input.stick));
    printf(__FUNCTION__ " : JoystickNumHats   : %i\n",
           SDL_JoystickNumHats(input.stick));
  } else if (xbox360) {
    // find and open the first xbox 360 controller we see.
    SDL_InitSubSystem(SDL_INIT_JOYSTICK);
    for (int i = 0; (input.stick = SDL_JoystickOpen(i)) != NULL; ++i) {
      string name(SDL_JoystickName(input.stick));
      printf(__FUNCTION__ " : JoystickName '%s'\n", name.c_str());
      if (name.find("X") == 0)
        break; // got it
      SDL_JoystickClose(input.stick);
    }
  }

  // Set up the video mode, OpenGL state, shaders and shader parameters.
  if (config.fullscreen) {
    initGraphics(true, 0, 0, lifeform.empty());
  } else {
    initGraphics(false, config.width, config.height, lifeform.empty());
  }

  // Parse as many uniforms from glsl source as we can find.
  uniforms.parseFromGlsl(::render.shaderManager.glsl_source);

  // TODO: prune uniforms to just those reported active by shader compiler.
  cout << ::render.shaderManager.fractal.uniforms();

  // Bind as many uniforms as we can find a match for to camera.
  uniforms.link(&camera);

  // Load key frames, if any. Clobbers camera.
  LoadKeyFrames(fixedFov);

  // Load initial camera; sets all known, linked uniforms.
  camera.loadConfig(BaseFile);
  if (fixedFov) {
    camera.width = config.width;
    camera.height = config.height;
    camera.fov_x = config.fov_x;
    camera.fov_y = config.fov_y;
  }

  initTwBar(stereoMode);
  fpsCounter.init(FPS_FRAMES_TO_AVERAGE);

  // printf(__FUNCTION__ " : GL_EXTENSIONS: %s\n", glGetString(GL_EXTENSIONS));

  // Main loop.
  Controller ctl = CTL_CAM; // the default controller is camera rotation
  int consecutiveChanges = 0;

  int done = 0;
  int frameno = 0;

  vector<KeyFrame> splines;
  size_t splines_index = 0;

  double frame_time = 1 / config.fps;
  printf(__FUNCTION__ " : frame time %g\n", frame_time);
  double render_time = 0;
  double render_start = 0;

  if (config.loop && !splines.empty()) {
    render_start = fpsCounter.now();
  }

  double dist_along_spline = 0;
  size_t keyframe = keyframes.size();

  bool ignoreNextMouseUp = false;
  bool dragging = false;
  struct {
    GLSL::dvec3 center;
    double dist;
    float step_x;
    float step_y;
  } DragCtx;
  bool pausing = false;
  bool stepping = false;

  bool mixedInOculus = false;
  bool multiPass = false;

  GLSL::vec3 effects_zoom(.5, .5, 1.); // centered and no zoom
  int zoomMouseX = 0;
  int zoomMouseY = 0;

  if (rendering) {
    // Rendering a sequence to disk. Spline the keyframes now.
    CatmullRom(keyframes, &splines, config.loop);
  }

  double last_de = 10.0;
  KeyFrame *next_camera = &camera;

  // Check availability of DE; setup func ptr.
  if (!de_func_name.empty()) {
    if (DE64_funcs.find(de_func_name) != DE64_funcs.end()) {
      // First pick is double version.
      de_func_64 = DE64_funcs[de_func_name];
    } else if (DE_funcs.find(de_func_name) != DE_funcs.end()) {
      // Next pick float version.
      de_func = DE_funcs[de_func_name];
    } else {
      printf(__FUNCTION__ " : unknown DE %s\n", de_func_name.c_str());
      de_func_name.clear();
    }
  }

  float view_q[4] = {0, 0, 0, 1};

  uint32_t poll_ticks = SDL_GetTicks();

  while (!done) {
    int ctlXChanged = 0, ctlYChanged = 0;

    if (!pausing || stepping) {

      if (next_camera != &camera)
        camera = *next_camera;

      // Set up current camera from spline.
      // Splined keyframes playback logic. Messy.
      if (!splines.empty()) {
        if (rendering && splines_index >= splines.size())
          break; // done

        // Loop if asked to.
        if (config.loop && splines_index >= splines.size()) {
          splines_index = 0;
          render_time = 0;
          render_start = fpsCounter.now();
        }

        // Figure out whether to draw a splined frame or skip to next.
        if (splines_index < splines.size()) {
          if (!config.loop && splines_index == splines.size() - 1) {
            camera = keyframes[keyframes.size() - 1]; // End at last position.
          } else {
            camera = splines[splines_index];
          }
          if (splines_index > 0) {
            dist_along_spline += camera.distanceTo(splines[splines_index - 1]);
          }
          size_t prev_splines_index = splines_index;
          ++splines_index;
          if (useTime) {
            if (rendering) {
              // Rendering sequece, use desired fps timing.
              if (camera.time < render_time)
                continue;
              render_time += frame_time;
            } else {
              // Previewing. Use real time (low framerate == jumpy preview!).
              double n = fpsCounter.now();
              if (n > render_start + camera.time)
                continue; // late, skip frame.
              double w = (render_start + camera.time) - n;
              if (w >= frame_time) { // early, redraw frame.
                splines_index = prev_splines_index;
              }
            }
          } else {
            if (dist_along_spline < (configSpeed ? config.speed : camera.speed))
              continue;
          }
          dist_along_spline -= (configSpeed ? config.speed : camera.speed);
        } else {
          splines.clear();
        }
      } else {
        camera.time = fpsCounter.now();
      }

      // We want current camera to always be in sync with some fields from
      // global config.
      camera.width = config.width;
      camera.height = config.height;
      camera.ipd = config.ipd;

      next_camera = &camera;

      if (!rendering) {
#if defined(_WIN32)
        // When not rendering a sequence,
        // now mix in orientation (and translation.. where's my DK2 Oculus?)
        if (stereoMode == ST_OCULUS) {
          if (GetOculusQuat(view_q)) {
            if (input.grabbed == true) {
              camera.mixSensorOrientation(view_q);
              mixedInOculus = true;
            }
          }
        }
#endif

        if (!config.disable_de) {
          // Try get a DE for current position.
          // We use DE as speed and interpupillairy distance (aka head size).
          double de = 0.0;
          if (de_func || de_func_64) {
            // Get a distance estimate, for navigation and eye separation.
            // Setup just the vars needed for DE. For now, iters and par[0..20]
            // TODO: make de() method of camera?
            GLSL::iters = camera.iters;
            GLSL::julia = camera.julia;
            for (int i = 0; i < 20; ++i) {
              GLSL::par[i] = GLSL::vec3(camera.par[i]);
            }

            GLSL::dvec3 pos(camera.pos());
            de = de_func_64 ? GLSL::abs(de_func_64(pos))
                            : GLSL::abs(de_func(pos));

          } else if (::render.shaderManager.de_shader.ok()) {
            // We did not find a CPU based DE; try the GPU based one.
            glDisable(GL_DEPTH_TEST);
            glBindFramebuffer(GL_FRAMEBUFFER, render.de_fbo);

            // Read previous' round computation, stored as pixel.
            GLSL::vec3 col = getPixelColor(0, 0);
            de = fabs(col.x);

            // Compute DE using shader.
            // We'll read the result next round so no costly glFinish needed.
            // TODO: use for auto-focus using center-weighted samples?
            GLuint program = ::render.shaderManager.de_shader.program();
            glUseProgram(program);
            camera.render(ST_COMPUTE_DE_ONLY, VQ_FRONT, program,
                          ::render.polarity);

            glUseProgram(0);
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
          }

          if (de != 0.0 && de != last_de) {
            printf("de=%12.12e\n", de);
            camera.speed = de / 10.0;
            last_de = de;
          }
        } // !config.disable_de
      }   // !rendering

      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_ALWAYS); // we're writing Z every pixel
      // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      GLuint program = ::render.shaderManager.fractal.program();
      glUseProgram(program); // the ::render.shaderManager.fractal shader

      // Figure out whether to render direct or via fbo.
      // A backbuffer (e.g. previous frame), or ::render.shaderManager.dof or
      // ::render.shaderManager.fxaa requires rendering to fbo.
      multiPass = (config.enable_dof &&
                   ((::render.shaderManager.dof.ok() && camera.enable_dof &&
                     camera.aperture != 0) ||
                    (::render.shaderManager.fxaa.ok() && camera.fxaa))) ||
                  config.backbuffer;

      // Set up input texture to ::render.shaderManager.fractal shader.
      if (render.background_texture) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, render.background_texture);
        glUniform1i(glGetUniformLocation(program, "iChannel0"), 0);
      }

      // Set up backbuffer as input.
      if (config.backbuffer) {
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, render.mainTex[(frameno & 1) ^ 1]);
        glUniform1i(glGetUniformLocation(program, "iBackbuffer"), 1);
        glUniform1i(glGetUniformLocation(program, "iBackbufferCount"),
                    camera.iBackbufferCount);
      }

      glUniform1i(glGetUniformLocation(program, "frameno"), frameno);

      if (multiPass) {
        glBindFramebuffer(GL_FRAMEBUFFER, render.mainFbo[frameno & 1]);
      }

      camera.speed *= speed_factor;

      if (rendercubes) {
        glViewport(0, 0, config.width, config.height);

        for (int vq = VQ_FRONT; vq < VQ_DONE; ++vq) {
          // Hack: render to render.fxaaFbo
          glBindFramebuffer(GL_FRAMEBUFFER, render.mainFbo[frameno & 1]);

          camera.render(stereoMode, vq, program, ::render.polarity);

          char filename[256];
          sprintf(filename, "cube-%05d%c.tga", frameno, VQ_LETTER[vq]);
          saveScreenshot(filename);
        }

        glViewport(0, 0, window.width(), window.height());
      } else {
        if (stereoMode == ST_DOME) {
          camera.render(stereoMode, VQ_UP, program, ::render.polarity);
        } else {
          // This is where the tflops go..
          camera.render(stereoMode, VQ_FRONT, program, ::render.polarity);
        }
      }

      camera.speed /= speed_factor;

      glUseProgram(0);

      glBindTexture(GL_TEXTURE_2D, 0);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);

    } // !pausing || stepping

    if (multiPass) {
      // We are post-processing of sorts.
      // Ortho projection; entire screen in regular pixel coordinates.
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(0, config.width, config.height, 0, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      GLuint currentFrame = render.mainTex[frameno & 1];

      if (!config.backbuffer) {
        glEnable(GL_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, currentFrame);
        // Mipmap the frame.
        // TODO: not needed by default. Costly?
        // glGenerateMipmap(GL_TEXTURE_2D);
      }

      if (::render.shaderManager.fxaa.ok() && config.enable_fxaa &&
          camera.fxaa) {
        // We have a ::render.shaderManager.fxaa shader.
        // Compute and point currentFrame(s) at output.
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, currentFrame);

        glBindFramebuffer(GL_FRAMEBUFFER, render.fxaaFbo);
        GLuint program = ::render.shaderManager.fxaa.program();
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "iTexture"), 0);
        glUniform1f(glGetUniformLocation(program, "xres"), config.width);
        glUniform1f(glGetUniformLocation(program, "yres"), config.height);
        drawScreen();

        glUseProgram(0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        currentFrame = render.fxaaTex; // Our output is input for next stage.
      }                                // ::render.shaderManager.fxaa

      if (::render.shaderManager.dof.ok() && config.enable_dof &&
          camera.enable_dof && camera.aperture != 0) {
        // We have a ::render.shaderManager.dof shader.
        // Compute iBlur0 and iBlur1
        for (int i = 0; i < NBLUR * 2; ++i) {

          // Setup input texture.
          int readLevel = 0;
          int writeLevel = 0;

          switch (i) {
          case 0: {
            // First input is current frame.
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, currentFrame);
          } break;
          default: {
            int odd = (i & 1);
            if (odd) {
              glActiveTexture(GL_TEXTURE0);
              glBindTexture(GL_TEXTURE_2D, render.scratchTex);
            } else {
              glActiveTexture(GL_TEXTURE0);
              glBindTexture(GL_TEXTURE_2D, currentFrame);
            }
            writeLevel = (i + 1) / 2 - 1;
            readLevel = (i / 2);
          } break;
          }

          // setup output buffer
          glBindFramebuffer(
              GL_FRAMEBUFFER,
              ((i & 1) == 0)
                  ? render.scratchFbo            // scratch as output
                  : render.blurFbo[writeLevel]); // target blur as output

          glActiveTexture(GL_TEXTURE2);
          glBindTexture(GL_TEXTURE_2D,
                        render.mainDepth[frameno & 1]); // current depth

          GLuint dof_program = ::render.shaderManager.dof.program();
          glUseProgram(dof_program);

          glUniform1i(glGetUniformLocation(dof_program, "iTexture"), 0);
          glUniform1i(glGetUniformLocation(dof_program, "iDepth"), 2);

          glUniform1i(glGetUniformLocation(dof_program, "iXorY"), (i & 1));

          glUniform1i(glGetUniformLocation(dof_program, "iBlurLevel"),
                      readLevel);

          glUniform1f(glGetUniformLocation(dof_program, "xres"), config.width);
          glUniform1f(glGetUniformLocation(dof_program, "yres"), config.height);

          drawScreen();

          glUseProgram(0); // no program
          // glDisable(GL_TEXTURE_2D);  // no texture
          glBindFramebuffer(GL_FRAMEBUFFER, 0); // default framebuffer
        }
      } // dof

      // Combine input(s) into final frame.
      GLuint final_program = ::render.shaderManager.effects.program();
      glUseProgram(final_program);

      glEnable(GL_TEXTURE_2D);

      // Pass main render as texture 0
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, currentFrame);

      if (!lifeform.empty()) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      }

      glUniform1i(glGetUniformLocation(final_program, "iTexture"), 0);

      // Pass main depth as texture 1
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, render.mainDepth[frameno & 1]);
      glUniform1i(glGetUniformLocation(final_program, "iDepth"), 1);

      if (::render.shaderManager.dof.ok() && camera.enable_dof &&
          camera.aperture != 0) {
        // Pass blur textures as 2 and 3, if we computed them.
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, render.blurTex[0]);
        glUniform1i(glGetUniformLocation(final_program, "iBlur0"), 2);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, render.blurTex[1]);
        glUniform1i(glGetUniformLocation(final_program, "iBlur1"), 3);
        glUniform1i(glGetUniformLocation(final_program, "enable_dof"), 1);
      } else {
        // No dof textures got computed.
        glUniform3f(glGetUniformLocation(final_program, "iZoom"),
                    effects_zoom.x, effects_zoom.y, effects_zoom.z);
        glUniform1i(glGetUniformLocation(final_program, "enable_dof"), 0);
      }

      // send tonemap params etc.
      glUniform1f(glGetUniformLocation(final_program, "exposure"),
                  config.exposure);
      glUniform1f(glGetUniformLocation(final_program, "maxBright"),
                  config.maxBright);
      glUniform1f(glGetUniformLocation(final_program, "gamma"), config.gamma);

      if (stereoMode == ST_OCULUS) {
        glUniform1f(glGetUniformLocation(final_program, "xres"), config.width);
        glUniform1f(glGetUniformLocation(final_program, "yres"), config.height);
        glUniform1f(glGetUniformLocation(final_program, "ipd"), config.ipd);
      }

      drawScreen();

      glUseProgram(0);
      // glDisable(GL_TEXTURE_2D);
    } // multiPass

    // Draw keyframe splined path, if we have 2+ keyframes and not rendering.
    if (!config.no_spline && stereoMode == ST_NONE && // TODO: show path in 3d
        keyframes.size() > 1 && splines.empty()) {

      // Wtf? need to deactivate each texture unit; glDisable is not enough.
      // Otherwise, texturing messes with color drawn and keyframe select fails.
      glEnable(GL_TEXTURE_2D);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, 0);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, 0);
      glActiveTexture(GL_TEXTURE2);
      glBindTexture(GL_TEXTURE_2D, 0);
      glActiveTexture(GL_TEXTURE3);
      glBindTexture(GL_TEXTURE_2D, 0);
      // glDisable(GL_TEXTURE_2D);

      glDepthFunc(GL_LESS);

      camera.activateGl();

      vector<KeyFrame> splines;
      CatmullRom(keyframes, &splines, config.loop);

      glColor4f(0.8f, 0.8f, 0.8f, 1.0f);
      glLineWidth(1);
      glEnable(GL_LINE_SMOOTH);

      glBegin(config.loop ? GL_LINE_LOOP : GL_LINE_STRIP); // right eye
      for (size_t i = 0; i < splines.size(); ++i) {
        splines[i].move(splines[i].speed, 0, 0);
        // Translate points, eye is at origin.
        glVertex3d(splines[i].pos()[0] - camera.pos()[0],
                   splines[i].pos()[1] - camera.pos()[1],
                   splines[i].pos()[2] - camera.pos()[2]);

        splines[i].move(-2 * splines[i].speed, 0, 0);
      }
      glEnd();
      glBegin(config.loop ? GL_LINE_LOOP : GL_LINE_STRIP); // left eye
      for (size_t i = 0; i < splines.size(); ++i) {
        // Translate points, eye is at origin.
        glVertex3d(splines[i].pos()[0] - camera.pos()[0],
                   splines[i].pos()[1] - camera.pos()[1],
                   splines[i].pos()[2] - camera.pos()[2]);
      }
      glEnd();

#if 0
      // draw light 1 path
      glColor4f(.1,.9,.9,1);
      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // light
      for (size_t i = 0; i < splines.size(); ++i) {
        // glVertex3dv(splines[i].par[19]);
        glVertex3d(splines[i].par[19][0] - camera.pos()[0],
                   splines[i].par[19][1] - camera.pos()[1],
                   splines[i].par[19][2] - camera.pos()[2]);
      }
      glEnd();
#endif

      glLineWidth(13);
      for (size_t i = 0; i < keyframes.size(); ++i) {
        // Encode keyframe # in blue channel.
        glColor4f(0.0f, 0.0f, 1.0f - (i / 256.0), 1.0f);
        glBegin(GL_LINES);
        KeyFrame tmp = keyframes[i];
        tmp.move(tmp.speed, 0, 0);
        // Translate points, eye is at origin.
        glVertex3d(tmp.pos()[0] - camera.pos()[0],
                   tmp.pos()[1] - camera.pos()[1],
                   tmp.pos()[2] - camera.pos()[2]);
        tmp = keyframes[i];
        tmp.move(-tmp.speed, 0, 0);
        glVertex3d(tmp.pos()[0] - camera.pos()[0],
                   tmp.pos()[1] - camera.pos()[1],
                   tmp.pos()[2] - camera.pos()[2]);
        glEnd();
      }
    }

    glDisable(GL_DEPTH_TEST);

    CHECK_ERROR;

    // Draw AntTweakBar
    if (input.grabbed == false && !rendering) {
      // glBindFramebuffer(GL_FRAMEBUFFER, 0);
      // glDisable(GL_TEXTURE_2D);
      // glDisable(GL_LINE_SMOOTH);
      TwDraw();
    }

    CHECK_ERROR;

    SDL_GL_SwapWindow(window.window());

    if (!pausing || stepping) {
      frameno++;
      camera.iBackbufferCount++;
      stepping = false;
    }

    if (rendering && !rendercubes && !splines.empty()) {
      // If we're playing a sequence back, save every frame to disk.
      char filename[256];
      sprintf(filename, "frame-%05d.tga", frameno);

      glBindFramebuffer(GL_FRAMEBUFFER, GL_FRONT);
      saveScreenshot(filename);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    fpsCounter.update();

    // Show position and fps in the caption.
    char caption[2048], controllerStr[256];
    sprintf(caption, "%s %.2ffps %5lu [%.3lf %.3lf %.3lf] %dms %s",
            printController(controllerStr, ctl), fpsCounter.get(),
            (unsigned long)splines_index, camera.pos()[0], camera.pos()[1],
            camera.pos()[2], fpsCounter.getLastFrameDuration(),
            pausing ? "(paused)" : "");
    SDL_SetWindowTitle(window.window(), caption);

    // Process UI events.
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_JOYBUTTONDOWN) {
        // Hack to reach "case SDLK_*" code:
        if (event.jbutton.button & 1)
          event.key.keysym.sym = SDLK_TAB;
        else
          event.key.keysym.sym = SDLK_BACKSPACE;
        event.key.keysym.mod = 0;
        event.type = SDL_KEYDOWN;
      }

#if defined(_WIN32)
      // Save current cursor in case Tw*() mucks with it.
      HCURSOR old_cur = ::GetCursor();
#endif
      if (input.grabbed == true ||
          !TwEventSDL(&event, SDL_MAJOR_VERSION, SDL_MINOR_VERSION)) {
#if defined(_WIN32)
        // Restore saved cursor after Tw*() call.
        if (::GetCursor() != old_cur) {
          ::SetCursor(old_cur);
        }
#endif
        switch (event.type) {
        case SDL_WINDOWEVENT: {
          if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
            if (initGraphics(false, event.window.data1, event.window.data2,
                             lifeform.empty())) {
              frameno = 0;
              stepping = true;
              initTwBar(stereoMode);

              config.fov_x = 0; // go for square pixels..
              config.sanitizeParameters();
            }
          }
        } break;

        case SDL_QUIT:
          done |= 1;
          break;

        case SDL_MOUSEBUTTONDOWN: {
          if (input.grabbed == false) {
            unsigned int bgr = getBGRpixel(event.button.x, event.button.y);
            if ((bgr & 0xffff) == 0) {
              // No red or green at all : probably a keyframe marker (fragile).
              size_t kf = 255 - (bgr >> 16);
              if (kf < keyframes.size()) {
                printf("selected keyframe %lu\n", (unsigned long)kf);
                keyframe = kf;
                dragging = true;
                ignoreNextMouseUp = true;
                XXX_SDL_SetCursor(input.hand);

                // Compute center point of dragging plane,
                // the plane though keyframe and orthogonal to view dir.
                double raw_dist = camera.distanceTo(keyframes[keyframe]);
                float dx = -event.button.x + config.width / 2;  // [-w/2 .. w/2]
                float dy = -config.height / 2 + event.button.y; // [-h/2 .. h/2]
                float ax = (camera.fov_x / (float)config.width) * dx;
                float ay = (camera.fov_y / (float)config.height) * dy;
                // Straight ahead dist to keyframe plane
                double dist =
                    raw_dist * cos(GLSL::radians(ax)) * cos(GLSL::radians(ay));

                GLSL::dvec3 O(camera.pos()), dir(camera.ahead());
                O += dir * dist; // Center on keyframe plane, parallel to screen

                // Compute step per pixel at distance 1
                float step_x =
                    tan(GLSL::radians(camera.fov_x / 2)) / (config.width / 2);
                float step_y =
                    tan(GLSL::radians(camera.fov_y / 2)) / (config.height / 2);

                // Populate DragCtx
                DragCtx.center = O;
                DragCtx.dist = dist;
                DragCtx.step_x = step_x;
                DragCtx.step_y = step_y;
              }
            }
          } else { // input.grabbed == true
            switch (event.button.button) {
            case 1: // left mouse
              input.grabbed = false;
              SDL_SetRelativeMouseMode(SDL_FALSE);
              ignoreNextMouseUp = true;
              break;
            case 4: // mouse wheel up, increase eye distance
              camera.speed *= 1.1;
              break;
            case 5: // mouse wheel down, decrease eye distance
              camera.speed *= 0.9;
              break;
            }
          }
        } break;

        case SDL_MOUSEMOTION: {
          if (effects_zoom.z == 1.0) {
            // Only update zoom center if zoomed out.
            zoomMouseX = event.motion.x;
            zoomMouseY = event.motion.y;
          }
          if (input.grabbed == false) {
            if (dragging == false) {
              // Peek at framebuffer color for keyframe markers.
              unsigned int bgr = getBGRpixel(event.motion.x, event.motion.y);
              // printf("bgr = %06x\n", bgr);
              if ((bgr & 0xffff) == 0) {
                // No red or green at all : probably a keyframe marker.
                size_t kf = 255 - (bgr >> 16);
                if (kf < keyframes.size()) {
                  printf("keyframe %lu\n", (unsigned long)kf);
                  if (SDL_GetCursor() == input.arrow) {
                    XXX_SDL_SetCursor(input.hand);
                  }
                } else {
                  if (SDL_GetCursor() == input.hand) {
                    XXX_SDL_SetCursor(input.arrow);
                  }
                }
              } else {
                if (SDL_GetCursor() == input.hand) {
                  XXX_SDL_SetCursor(input.arrow);
                }
              }
            } else { // dragging == true
              if (SDL_GetCursor() != input.hand) {
                XXX_SDL_SetCursor(input.hand);
              }
              // Drag currently selected keyframe around.
              if (keyframe < keyframes.size()) {
                GLSL::dvec3 v(DragCtx.center);
                GLSL::dvec3 up(camera.up()), right(camera.right());
                float dx = -event.motion.x + config.width / 2;  // [-w/2 .. w/2]
                float dy = -config.height / 2 + event.motion.y; // [-h/2 .. h/2]
                                                                //
                v += right * DragCtx.step_x * -dx * DragCtx.dist;
                v += up * DragCtx.step_y * -dy * DragCtx.dist;

                keyframes[keyframe].set(v.x, v.y, v.z);
              }
            }
          }
        } break;

        case SDL_MOUSEBUTTONUP: {
          if (ignoreNextMouseUp == false && input.grabbed == false) {
            input.grabbed = true;
            if (lifeform.empty()) {
              XXX_SDL_SetCursor(input.arrow);
              SDL_SetRelativeMouseMode(SDL_TRUE);
            } else {
              XXX_SDL_SetCursor(input.crosshair);
            }
          }
          ignoreNextMouseUp = false;
          dragging = false;
        } break;

        case SDL_MOUSEWHEEL: {
          if (!lifeform.empty()) {
            effects_zoom.x = (float)zoomMouseX / config.width;
            effects_zoom.y =
                (float)(config.height - zoomMouseY) / config.height;
            effects_zoom.z =
                GLSL::clamp(effects_zoom.z + event.wheel.y, 1.f, 16.f);
          }
        } break;

        case SDL_KEYDOWN: {
          bool hasAlt = event.key.keysym.mod & (KMOD_LALT | KMOD_RALT);
          bool hasCtrl = event.key.keysym.mod & (KMOD_LCTRL | KMOD_RCTRL);
          bool hasShift = event.key.keysym.mod & (KMOD_LSHIFT | KMOD_RSHIFT);

          switch (event.key.keysym.sym) {
          case SDLK_ESCAPE: {
            if (input.grabbed == true) {
              input.grabbed = false;
              if (!lifeform.empty()) {
                SDL_SetRelativeMouseMode(SDL_FALSE);
              }
            } else {
              done |= 1;
            }
          } break;

          // Switch fullscreen mode (drops the whole OpenGL context in Windows).
          case SDLK_RETURN:
          case SDLK_KP_ENTER: {
            initGraphics(true, 0, 0, lifeform.empty());
            frameno = 0;
            stepping = true;
            initTwBar(stereoMode);
          } break;

          // Switch L/R eye polarity,
          // or print screenshot (w/ ctrl).
          case SDLK_p: {
            if (!hasCtrl) {
              ::render.polarity *= -1;
            } else {
              // Save config and screenshot (filename = current time).
              time_t t = time(0);
              struct tm *ptm = localtime(&t);
              char filename[256];
              strftime(filename, 256, "%Y%m%d_%H%M%S.cfg", ptm);
              camera.saveConfig(filename);
              strftime(filename, 256, "%Y%m%d_%H%M%S.tga", ptm);
              glBindFramebuffer(GL_FRAMEBUFFER, GL_FRONT);
              saveScreenshot(filename);
              glBindFramebuffer(GL_FRAMEBUFFER, 0);
            }
          } break;

          // Splined path control.
          case SDLK_HOME: { // Start spline playback from start.
            if (!keyframes.empty()) {
              CatmullRom(keyframes, &splines, config.loop);
              splines_index = 0;
              dist_along_spline = 0;
              render_time = 0;
              render_start = fpsCounter.now();
              keyframe = 0;
            } else
              next_camera = &config;
          } break;

          case SDLK_END: { // Stop spline playback.
            splines.clear();
            if (!keyframes.empty()) {
              keyframe = keyframes.size() - 1;
              next_camera = &keyframes[keyframe];
            } else
              next_camera = &config;
          } break;

          case SDLK_DELETE: {
            // Drop last keyframe, reset camera to previous keyframe.
            splines.clear();
            if (!keyframes.empty()) {
              if (hasCtrl) {
                // delete current keyframe
                if (keyframe < keyframes.size()) {
                  keyframes.erase(keyframes.begin() + keyframe);
                  if (keyframe >= keyframes.size())
                    keyframe = keyframes.size() - 1;
                }
              } else {
                // delete last keyframe
                keyframes.pop_back();
                keyframe = keyframes.size() - 1;

                if (keyframe < keyframes.size())
                  next_camera = &keyframes[keyframe];
              }
            }
          } break;

          case SDLK_SPACE:
            if (pausing) {
              stepping = true;
              break;
            }
          case SDLK_INSERT: { // Add keyframe.
            splines.clear();
            size_t index = keyframes.size();
            if (hasCtrl) {
              // Replace currently selected keyframe.
              index = keyframe;
            }
            if (!hasCtrl) {
              // Need an estimate for delta_time for this new frame.
              suggestDeltaTime(camera, keyframes);
              keyframes.push_back(camera); // Add keyframe at end.
            } else {
              keyframes[index] = camera; // Overwrite current keyframe.
            }
            char filename[256];
            sprintf(filename, "%s-%lu.cfg", kKEYFRAME, (unsigned long)index);
            camera.saveConfig(filename);
          } break;

          case SDLK_PAGEUP: { // Jump keyframe back in spline playback.
            if (splines.empty()) {
              CatmullRom(keyframes, &splines, config.loop);
              splines_index = splines.size();
            }
            if (splines_index < NSUBFRAMES)
              splines_index = 0;
            else
              splines_index -= NSUBFRAMES;
          } break;

          case SDLK_PAGEDOWN: { // Jump splined playback keyframe ahead.
            splines_index += NSUBFRAMES;
            if (splines_index >= splines.size()) {
              splines_index = splines.size() - 1;
            }
          } break;

          case SDLK_TAB:
            if (!hasShift) {
              { // Step thru keyframes.
                if (!splines.empty()) {
                  // find most recently played keyframe
                  int nkeys = 0;
                  for (size_t i = 0; i < splines_index; ++i)
                    if (splines[i].isKey())
                      ++nkeys;
                  keyframe = nkeys;
                }
                if (!hasCtrl)
                  ++keyframe;
                if (keyframe >= keyframes.size())
                  keyframe = 0;

                if (keyframe < keyframes.size() && splines.empty()) {
                  // Don't jump camera ahead if we were playing, just stop in
                  // place.
                  next_camera = &keyframes[keyframe];
                }

                if (keyframe < keyframes.size()) {
                  printf("at keyframe %lu, speed %.8e, delta_time %f\n",
                         (unsigned long)keyframe, keyframes[keyframe].speed,
                         keyframes[keyframe].delta_time);
                }

                if (keyframes.empty())
                  next_camera = &config; // back to start

                if (hasCtrl) {
                  // Start playing: spline and start at keyframe.
                  if (!keyframes.empty()) {
                    CatmullRom(keyframes, &splines, config.loop);
                    size_t nkeys = keyframe;
                    for (splines_index = 0; splines_index < splines.size();
                         ++splines_index) {
                      if (splines[splines_index].isKey())
                        if (nkeys-- == 0)
                          break;
                    }
                    render_time = splines[splines_index].time;
                    render_start = fpsCounter.now() - render_time;
                  }
                } else {
                  splines.clear();
                }
              }
              break;
            } // shift-tab == backspace, fall through

          case SDLK_BACKSPACE: {
            --keyframe;
            if (keyframe >= keyframes.size())
              keyframe = keyframes.size() - 1;
            if (keyframe < keyframes.size()) {
              next_camera = &keyframes[keyframe];
              printf("at keyframe %lu, speed %.8e, delta_time %f\n",
                     (unsigned long)keyframe, keyframes[keyframe].speed,
                     keyframes[keyframe].delta_time);
            } else
              next_camera = &config;
            splines.clear();
          } break;

          // Resolve controller value changes that happened during rendering.
          case SDLK_LEFT:
            ctlXChanged = 1;
            updateControllerX(ctl, -(consecutiveChanges = 1), hasAlt);
            break;
          case SDLK_RIGHT:
            ctlXChanged = 1;
            updateControllerX(ctl, (consecutiveChanges = 1), hasAlt);
            break;
          case SDLK_DOWN:
            ctlYChanged = 1;
            updateControllerY(ctl, -(consecutiveChanges = 1), hasAlt);
            break;
          case SDLK_UP:
            ctlYChanged = 1;
            updateControllerY(ctl, (consecutiveChanges = 1), hasAlt);
            break;

          // Currently selected keyframe manouvering in screenspace.
          // 4-6 left-right, 8-2 up-down, 9-1 further-closer, 7-3 speed at
          // frame.
          case SDLK_KP_4: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.right(),
                                               -.5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_6: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.right(),
                                               .5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_8: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.up(),
                                               .5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_2: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.up(),
                                               -.5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_9: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.ahead(),
                                               .5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_1: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].moveAbsolute(camera.ahead(),
                                               -.5 * keyframes[keyframe].speed);
            }
          } break;
          case SDLK_KP_3: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].speed *= 1.1;
            }
          } break;
          case SDLK_KP_7: {
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].speed *= .9;
            }
          } break;

          case SDLK_SCROLLLOCK:
          case SDLK_PAUSE: {
            pausing = !pausing;
          } break;

          // See whether the active controller has changed.
          default: {
            Controller oldCtl = ctl;
            changeController(event.key.keysym.sym, &ctl);
            if (ctl != oldCtl) {
              consecutiveChanges = 0;
            }
          } break;
          }
        } break;
        }
      }
    }

    if (done)
      break;

    // Poll keyboard and mouse state for level-based inputs.
    // At around 1000/16 per second polling rate (or slower if FPS is slower..).
    if (SDL_GetTicks() - poll_ticks >= 16) {
      poll_ticks += 16;

      const Uint8 *keystate = SDL_GetKeyboardState(0);
      int mouse_dx, mouse_dy;
      Sint16 joystick_x = 0, joystick_y = 0, joystick_z = 0, joystick_r = 0,
             joystick_lt = -32768, joystick_rt = -32768;
      Uint8 joystick_hat = 0;
      Uint8 mouse_buttons = SDL_GetRelativeMouseState(&mouse_dx, &mouse_dy);
      int mouse_button_left = mouse_buttons & SDL_BUTTON(SDL_BUTTON_LEFT);
      int mouse_button_right = mouse_buttons & SDL_BUTTON(SDL_BUTTON_RIGHT);

      bool hasAlt = keystate[SDL_SCANCODE_RALT] || keystate[SDL_SCANCODE_LALT];
      bool hasCtrl =
          keystate[SDL_SCANCODE_RCTRL] || keystate[SDL_SCANCODE_LCTRL];
      bool hasShift =
          keystate[SDL_SCANCODE_RSHIFT] || keystate[SDL_SCANCODE_LSHIFT];
      (void)hasAlt;
      (void)hasShift;
      (void)hasCtrl;

      // Continue after calling SDL_GetRelativeMouseState() so view direction
      // does not jump after closing AntTweakBar.
      if (input.grabbed == false) {
        if (mixedInOculus) {
          camera.unmixSensorOrientation(view_q);
        }
        continue;
      }

      if (input.stick) {
        SDL_JoystickUpdate();
        joystick_x = SDL_JoystickGetAxis(input.stick, 2);
        joystick_y = SDL_JoystickGetAxis(input.stick, 3);
        joystick_z = SDL_JoystickGetAxis(input.stick, 1);
        joystick_r = SDL_JoystickGetAxis(input.stick, 0);
        joystick_lt = SDL_JoystickGetAxis(input.stick, 4);
        joystick_rt = SDL_JoystickGetAxis(input.stick, 5);
        joystick_hat = SDL_JoystickGetHat(input.stick, 0);
        if (abs(joystick_x) < 5000)
          joystick_x = 0;
        if (abs(joystick_y) < 5000)
          joystick_y = 0;
        if (abs(joystick_z) < 5000)
          joystick_z = 0;
        if (abs(joystick_r) < 10000)
          joystick_r = 0;
      }

      (void)mouse_buttons;
      (void)mouse_button_left;
      (void)mouse_button_right;
      (void)joystick_lt;
      (void)joystick_rt;

      if (keystate[SDL_SCANCODE_W])
        camera.move(0, 0, camera.speed); // forward
      if (keystate[SDL_SCANCODE_S])
        camera.move(0, 0, -camera.speed); // back
#if !defined(MOVE_W_TRIGGERS)
      if (joystick_z != 0)
        camera.move(0, 0, camera.speed * -joystick_z / 20000.0);
#else
      if (joystick_lt != -32768)
        camera.move(0, 0, camera.speed * (joystick_lt + 32768) / 20000.0);
      if (joystick_rt != -32768)
        camera.move(0, 0, -camera.speed * (joystick_rt + 32768) / 20000.0);
#endif

      if (keystate[SDL_SCANCODE_A] || (joystick_hat & SDL_HAT_LEFT))
        camera.move(-camera.speed, 0, 0); // left
      if (keystate[SDL_SCANCODE_D] || (joystick_hat & SDL_HAT_RIGHT))
        camera.move(camera.speed, 0, 0); // right
      if (joystick_hat & SDL_HAT_DOWN)
        camera.move(0, -camera.speed, 0); // down
      if (joystick_hat & SDL_HAT_UP)
        camera.move(0, camera.speed, 0); // up

      // Mouse look.
      if (input.grabbed == true && (mouse_dx != 0 || mouse_dy != 0)) {
        m_rotateX2(camera.mouse_rot_speed * mouse_dx * camera.fov_x / 90.0);
        m_rotateY2(camera.mouse_rot_speed * mouse_dy * camera.fov_y / 75.0);
      }
      if (keystate[SDL_SCANCODE_Q])
        m_rotateZ2(camera.keyb_rot_speed);
      if (keystate[SDL_SCANCODE_E])
        m_rotateZ2(-camera.keyb_rot_speed);

      // input.stick look.
      if (joystick_x != 0 || joystick_y != 0) {
        m_rotateX2(camera.mouse_rot_speed * joystick_x * camera.fov_x / 90.0 /
                   20000.0);
        m_rotateY2(camera.mouse_rot_speed * -joystick_y * camera.fov_y / 75.0 /
                   10000.0);
      }
      if (joystick_r != 0)
        m_rotateZ2(camera.keyb_rot_speed * -joystick_r / 100000.0);

      if (keystate[SDL_SCANCODE_Z]) {
        if (camera.speed > 0.000001)
          camera.speed -= camera.speed / 10;
        printf("speed %.8e\n", camera.speed);
      }
      if (keystate[SDL_SCANCODE_C]) {
        if (camera.speed < 1.0)
          camera.speed += camera.speed / 10;
        printf("speed %.8e\n", camera.speed);
      }

      // Change the value of the active controller.
      if (!ctlXChanged) {
        if (keystate[SDL_SCANCODE_LEFT]) {
          ctlXChanged = 1;
          updateControllerX(ctl, -++consecutiveChanges, hasAlt);
        }
        if (keystate[SDL_SCANCODE_RIGHT]) {
          ctlXChanged = 1;
          updateControllerX(ctl, ++consecutiveChanges, hasAlt);
        }
      }
      if (!ctlYChanged) {
        if (keystate[SDL_SCANCODE_DOWN]) {
          ctlYChanged = 1;
          updateControllerY(ctl, -++consecutiveChanges, hasAlt);
        }
        if (keystate[SDL_SCANCODE_UP]) {
          ctlYChanged = 1;
          updateControllerY(ctl, ++consecutiveChanges, hasAlt);
        }
      }

#if defined(HYDRA)
      // Sixense Hydra
      if (sixenseGetAllNewestData(&ssdata) == SIXENSE_SUCCESS) {
        //  printf("%f %f %f %f\n",
        //  ssdata.controllers[0].rot_quat[0],ssdata.controllers[0].rot_quat[1],ssdata.controllers[0].rot_quat[2],ssdata.controllers[0].rot_quat[3]);

        int clbuttons = ssdata.controllers[1].buttons;
        int crbuttons = ssdata.controllers[0].buttons;

#if 0
  camera.move(0, 0, camera.speed *   ssdata.controllers[0].joystick_y);
  m_rotateX2(camera.keyb_rot_speed *.1 * ssdata.controllers[1].joystick_x);
  m_rotateY2(camera.keyb_rot_speed *.1 * ssdata.controllers[1].joystick_y);
  m_rotateZ2(camera.keyb_rot_speed *.1 * -ssdata.controllers[0].joystick_x);
#endif

        //  printf("%08x, %f\n", ssdata.controllers[0].buttons,
        //  ssdata.controllers[0].trigger); printf("%+7.7f %+7.7f %+7.7f, ",
        //  ssdata.controllers[0].pos[0],ssdata.controllers[0].pos[1],ssdata.controllers[0].pos[2]);
        //  printf("%+7.7f %+7.7f %+7.7f\n",
        //  ssdata.controllers[1].pos[0],ssdata.controllers[1].pos[1],ssdata.controllers[1].pos[2]);

        float dx = ssdata.controllers[0].pos[0] - ssdata.controllers[1].pos[0];
        float dy = ssdata.controllers[0].pos[1] - ssdata.controllers[1].pos[1];
        float dz = ssdata.controllers[0].pos[2] - ssdata.controllers[1].pos[2];
        float d = sqrt(dx * dx + dy * dy + dz * dz);

        // printf("%+7.7f\n", d);

        // spot in between two hands.
        dx =
            (ssdata.controllers[0].pos[0] + ssdata.controllers[1].pos[0]) / 2.0;
        dy =
            (ssdata.controllers[0].pos[1] + ssdata.controllers[1].pos[1]) / 2.0;
        dz =
            (ssdata.controllers[0].pos[2] + ssdata.controllers[1].pos[2]) / 2.0;

        // set neutral when any start button is pressed.
        bool calibrate = (lbuttons | rbuttons | clbuttons | crbuttons) &
                         SIXENSE_BUTTON_START;

        if (calibrate) {
          speed_base = d;
          neutral_x = dx;
          neutral_y = dy;
          neutral_z = dz;
        }

        // distance between controllers is eye separation / speed multiplier.
        speed_factor = d / speed_base;

        //  printf("triggers %+7.7f, %+7.7f\n", ssdata.controllers[0].trigger,
        //  ssdata.controllers[1].trigger); printf("buttons  %+7.7x, %+7.7x\n",
        //  clbuttons, crbuttons);

        // flip back&forth through keyframes w/ edge trigger of bumper button
        // presses.
        if ((clbuttons & SIXENSE_BUTTON_BUMPER) &&
            ((clbuttons ^ lbuttons) & SIXENSE_BUTTON_BUMPER)) {
          --keyframe;
          if (keyframe >= keyframes.size())
            keyframe = keyframes.size() - 1;
          if (keyframe < keyframes.size()) {
            next_camera = &keyframes[keyframe];
          } else {
            next_camera = &config;
          }
        }
        lbuttons = clbuttons;

        if ((crbuttons & SIXENSE_BUTTON_BUMPER) &&
            ((crbuttons ^ rbuttons) & SIXENSE_BUTTON_BUMPER)) {
          ++keyframe;
          if (keyframe >= keyframes.size())
            keyframe = 0;
          if (keyframe < keyframes.size()) {
            next_camera = &keyframes[keyframe];
          } else {
            next_camera = &config;
          }
        }
        rbuttons = crbuttons;

        dx = (neutral_x - dx) / 100.0;
        dy = (neutral_y - dy) / 100.0;
        dz = (neutral_z - dz) / 100.0;

        // printf("%+8.8lf %+8.8lf %+8.8lf\n", dx, dy, dz);

        dx *= (camera.speed * speed_factor);
        dy *= (camera.speed * speed_factor);
        dz *= (camera.speed * speed_factor);
        camera.move(-dx, -dy, dz);
      }
#endif // HYDRA

      if (!(ctlXChanged || ctlYChanged))
        consecutiveChanges = 0;

      // We might have changed view. Preserve changes, minus HMD orientation.
      if (mixedInOculus) {
        camera.unmixSensorOrientation(view_q);
      }
    }

    if (outputFilename != NULL) {
      glBindFramebuffer(GL_FRAMEBUFFER, GL_FRONT);
      saveScreenshot(outputFilename);
      break;
    }
    if (!rendering && rendercubes) {
      break;
    }
  }

  TwTerminate();

  camera.width = saveWidth;
  camera.height = saveHeight;
  camera.saveConfig("last.cfg",
                    &defines); // Save a config file on exit, just in case.

#if defined(_WIN32)
  if (stereoMode == ST_OCULUS) {
    ReleaseOculusSDK();
  }
#endif

#if defined(HYDRA)
  sixenseExit();
#endif

  window.reset();

  return 0;
}
