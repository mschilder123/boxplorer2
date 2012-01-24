#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>

#if !defined(_WIN32)

#include <float.h>
#include <unistd.h>
#define MKDIR(a) mkdir((a), 0755)

#define _strdup strdup
#define __FUNCTION__ "boxplorer"
#define MAX_PATH 256

#else  // _WIN32

#include <direct.h>
#define MKDIR(a) mkdir(a)

#pragma warning(disable: 4996) // unsafe function
#pragma warning(disable: 4244) // conversion loss
#pragma warning(disable: 4305) // truncation
#pragma warning(disable: 4800) // forcing value to bool

#pragma comment(lib, "SDL.lib")
#pragma comment(lib, "SDLmain.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "shell32.lib")

#endif  // _WIN32

#include <vector>
#include <string>
#include <map>

using namespace std;

#define NO_SDL_GLEXT
#include <SDL/SDL_opengl.h>
#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#include <SDL/SDL_main.h>

#include <AntTweakBar.h>

#include "shader_procs.h"
#include "default_shaders.h"

#include "TGA.h"

#define DEFAULT_CONFIG_FILE  "boxplorer.cfg"
#define VERTEX_SHADER_FILE   "vertex.glsl"
#define FRAGMENT_SHADER_FILE "fragment.glsl"
#define FRAME_VERTEX_SHADER_FILE   "frame_vertex.glsl"
#define FRAME_FRAGMENT_SHADER_FILE "frame_fragment.glsl"

#ifdef PI
  #undef PI
#endif
#define PI          3.14159265358979324
#define die(...)    ( fprintf(stderr, __VA_ARGS__), exit(-1), 1 )
#define lengthof(x) ( sizeof(x)/sizeof((x)[0]) )

#include "glsl.h"

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
    if (x64 != string::npos) name.erase(x64);
    DE64_funcs[name] = func;
  }
};
#define DECLARE_DE(a) DE_initializer _init##a(#a, &a);
#define DECLARE_COLORING(a)  // not interested in color funcs here

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
};

#define sign(a) GLSL::sign(a)

#define zNear 0.0001f
#define zFar  5.0f

#define FPS_FRAMES_TO_AVERAGE 20

const char* kKEYFRAME = "keyframe";

static const char *kHand[] = {
  "     XX                 ",
  "    X..X                ",
  "    X..X                ",
  "    X..X                ",
  "    X..XXXXXX           ",
  "    X..X..X..XXX        ",
  "XXX X..X..X..X..X       ",
  "X..XX..X..X..X..X       ",
  "X...X..X..X..X..X       ",
  " X..X..X..X..X..X       ",
  " X..X........X..X       ",
  " X..X...........X       ",
  " X..............X       ",
  " X.............X        ",
  " X.............X        ",
  "  X...........X         ",
  "   X.........X          ",
  "    X........X          ",
  "    X........X          ",
  "    XXXXXXXXXX          ",
};

static SDL_Cursor *init_system_cursor(const char *image[]) {
  int i = -1;
  Uint8 data[3*20];
  Uint8 mask[3*20];

  for ( int row=0; row<20; ++row ) {
    for ( int col=0; col<24; ++col ) {
      if ( col % 8 ) {
        data[i] <<= 1; mask[i] <<= 1;
      } else {
        ++i;
        data[i] = mask[i] = 0;
      }
      switch (image[row][col]) {
        case '.':
          data[i] |= 0x01; mask[i] |= 0x01;
          break;
        case 'X':
          mask[i] |= 0x01;
          break;
        case ' ':
          break;
      }
    }
  }
  return SDL_CreateCursor(data, mask, 24, 20, 5, 0);
}

// SDL cursors.
SDL_Cursor* arrow_cursor;
SDL_Cursor* hand_cursor;

// The shader program handle.
int program;
// Our SDL handle.
SDL_Surface* screen;
// Optional #defines for glsl compilation from .cfg file.
string defines;

// Pinhole camera modes.
enum StereoMode { ST_NONE=0, ST_OVERUNDER, ST_XEYED, ST_INTERLACED, ST_SIDEBYSIDE }
    stereoMode = ST_NONE;

// ogl framebuffer object
GLuint fbo;
// texture that frame got rendered to
GLuint texture;
// depth buffer attached to fbo
GLuint depthBuffer;
// the dof mipmapper program
int dof_program;
// texture holding background image
GLuint background_texture;

////////////////////////////////////////////////////////////////
// Helper functions

// Compute the dot product of two vectors.
double dot(const double x[3], const double y[3]) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// Normalize a vector. If it was zero, return 0.
int normalize(double x[3]) {
  double len = dot(x, x); if (len == 0) return 0;
  len = 1/sqrt(len); x[0] *= len; x[1] *= len; x[2] *= len;
  return 1;
}

// Allocate a char[] and read a text file into it. Return 0 on error.
char* readFile(char const* name) {
  FILE* f;
  int len;
  char* s = 0;

  // open file an get its length
  if (!(f = fopen(name, "r"))) goto readFileError1;
  fseek(f, 0, SEEK_END);
  len = ftell(f);

  // read the file in an allocated buffer
  if (!(s = (char*)malloc(len+1))) goto readFileError2;
  rewind(f);
  len = fread(s, 1, len, f);
  s[len] = '\0';

  readFileError2: fclose(f);
  readFileError1: return s;
}


////////////////////////////////////////////////////////////////
// FPS tracking.

int framesToAverage;
Uint32* frameDurations;
int frameDurationsIndex = 0;
Uint32 lastFrameTime;

float now() {
  return (float)SDL_GetTicks() / 1000.0;
}

// Initialize the FPS structure.
void initFPS(int framesToAverage_) {
  assert(framesToAverage_ > 1);
  framesToAverage = framesToAverage_;
  frameDurations = (Uint32*)malloc(sizeof(Uint32) * framesToAverage_);
  frameDurations[0] = 0;
  lastFrameTime = SDL_GetTicks();
}

// Update the FPS structure after drawing a frame.
void updateFPS(void) {
  Uint32 time = SDL_GetTicks();
  frameDurations[frameDurationsIndex++ % framesToAverage] =
      time - lastFrameTime;
  lastFrameTime = time;
}

// Return the duration of the last frame.
Uint32 getLastFrameDuration(void) {
  return frameDurations[(frameDurationsIndex+framesToAverage-1) %
      framesToAverage];
}

// Return the average FPS over the last X frames.
float getFPS(void) {
  if (frameDurationsIndex < framesToAverage) return 0;  // not enough data

  int i; Uint32 sum;
  for (i=sum=0; i<framesToAverage; i++) sum += frameDurations[i];
  return framesToAverage * 1000.f / sum;
}


////////////////////////////////////////////////////////////////
// Current logical state of the program.

// Configuration parameters.
// [type, variable name, config name, spline]
#define PROCESS_CONFIG_PARAMS \
  PROCESS(int, width, "width", false) \
  PROCESS(int, height, "height", false) \
  PROCESS(int, fullscreen, "fullscreen", false) \
  PROCESS(int, multisamples, "multisamples", false) \
  PROCESS(float, keyb_rot_speed, "keyb_rot_speed", false) \
  PROCESS(float, mouse_rot_speed, "mouse_rot_speed", false) \
  PROCESS(float, fov_x, "fov_x", true) \
  PROCESS(float, fov_y, "fov_y", true) \
  PROCESS(double, speed, "speed", true) \
  PROCESS(float, min_dist, "min_dist", true) \
  PROCESS(int, max_steps, "max_steps", true) \
  PROCESS(int, iters, "iters", true) \
  PROCESS(int, color_iters, "color_iters", true) \
  PROCESS(int, loop, "loop", false) \
  PROCESS(float, ao_eps, "ao_eps", true) \
  PROCESS(float, ao_strength, "ao_strength", true) \
  PROCESS(float, glow_strength, "glow_strength", true) \
  PROCESS(float, dist_to_color, "dist_to_color", true) \
  PROCESS(float, delta_time, "delta_time", false) \
  PROCESS(float, time, "time", true) \
  PROCESS(float, fps, "fps", false) \
  PROCESS(int, depth_size, "depth_size", false) \
  PROCESS(float, z_near, "z_near", true) \
  PROCESS(float, z_far, "z_far", true) \
  PROCESS(float, dof_scale, "dof_scale", true) \
  PROCESS(float, dof_offset, "dof_offset", true) \
  PROCESS(int, enable_dof, "enable_dof", false) \
  PROCESS(int, no_spline, "no_spline", false) \
  PROCESS(float, focus, "focus", true) \
  PROCESS(int, nrays, "nrays", true)

#define NUMPARS 20
char* parName[NUMPARS][3];

class KeyFrame {
  public:
   // View matrix.
   double v[16];

   // Declare fractal and other parameters.
   #define PROCESS(a,b,c,d) a b;
   PROCESS_CONFIG_PARAMS
   #undef PROCESS

   // Par[] parameter array.
   float par[NUMPARS][3];  // min(this, glsl) gets sent to shader.

   bool isKey_;  // Whether this frame is actually a KeyFrame.

   KeyFrame() { memset(this, 0, sizeof *this); }

   double distanceTo(const KeyFrame& other) const {
      double delta[3] = { v[12]-other.v[12],
                          v[13]-other.v[13],
                          v[14]-other.v[14] };
      return sqrt(dot(delta, delta));
   }
   double* right() { return &v[0]; }
   double* up() { return &v[4]; }
   double* ahead() { return &v[8]; }
   double* pos() { return &v[12]; }

   void setKey(bool key) { isKey_ = key; }
   bool isKey() const { return isKey_; }

   void orthogonalize() {
      if (!normalize(ahead())) { ahead()[0]=ahead()[1]=0; ahead()[2]=1; }
      double l = dot(ahead(), up());
      for (int i=0; i<3; i++) up()[i] -= l*ahead()[i];
      if (!normalize(up())) {  // Error? Make upDirection.z = 0.
         up()[2] = 0;
         if (fabs(ahead()[2]) == 1) {
            up()[0] = 0;
            up()[1] = 1;
         } else {
            up()[0] = -ahead()[1];
            up()[1] = ahead()[0];
            normalize(up());
         }
      }
      // Compute rightDirection as a cross product of upDirection and direction.
      for (int i=0; i<3; i++) {
         int j = (i+1)%3, k = (i+2)%3;
         right()[i] = up()[j]*ahead()[k] - up()[k]*ahead()[j];
      }
   }

   // Move camera in a direction relative to the view direction.
   // Behaves like `glTranslate`.
   void move(double x, double y, double z) {
      for (int i=0; i<3; i++) {
         pos()[i] += right()[i]*x + up()[i]*y + ahead()[i]*z;
      }
   }

   // Move camera in the normalized absolute direction `dir` by `len` units.
   void moveAbsolute(double* dir, double len) {
      for (int i=0; i<3; i++) {
         pos()[i] += len * dir[i];
      }
   }

   // Rotate the camera by `deg` degrees around a normalized axis.
   // Behaves like `glRotate` without normalizing the axis.
   void rotate(double deg, double x, double y, double z) {
     double s = sin(deg*PI/180), c = cos(deg*PI/180), t = 1-c;
     double r[3][3] = {
      { x*x*t +   c, x*y*t + z*s, x*z*t - y*s },
      { y*x*t - z*s, y*y*t +   c, y*z*t + x*s },
      { z*x*t + y*s, z*y*t - x*s, z*z*t +   c }
     };
     for (int i=0; i<3; i++) {
      double c[3];
      for (int j=0; j<3; j++) c[j] = v[i+j*4];
      for (int j=0; j<3; j++) v[i+j*4] = dot(c, r[j]);
     }
   }

   // Set the OpenGL modelview matrix to the camera matrix, for shader.
   void activate() const {
      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixd(v);
   }

   // Set the OpenGL modelview and projection for gl*() functions.
   void activateGl() {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      float fH = tan( fov_y * PI / 360.0f ) * zNear;
      float fW = tan( fov_x * PI / 360.0f ) * zNear;
      glFrustum(-fW, fW, -fH, fH, zNear, zFar);

      orthogonalize();
      double matrix[16] = {
         right()[0], up()[0], -ahead()[0], 0,
         right()[1], up()[1], -ahead()[1], 0,
         right()[2], up()[2], -ahead()[2], 0,
                  0,       0,           0, 1
      };
      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixd(matrix);
      glTranslated(-pos()[0], -pos()[1], -pos()[2]);
   }

   // Load configuration.
   bool loadConfig(char const* configFile, string* defines = NULL) {
     bool result = false;
     FILE* f;
     if ((f = fopen(configFile, "r")) != 0) {
      size_t i;
      char s[32768];  // max line length
      while (fscanf(f, " %s", s) == 1) {  // read word
        if (s[0] == 0 || s[0] == '#') continue;

        int v;

        // Parse #defines out of config.cfg to prepend to .glsl
        if (defines) {
          if (!strcmp(s, "d") || !strcmp(s, "c")) {
            string a(s);
            v = fscanf(f, " %s", s);
            if (v == 1) {
              string define = "#define " + a + " " + s + "\n";
              printf(__FUNCTION__ " : %s", define.c_str());
              defines->append(define);
              if (!a.compare("d")) {
                de_func_name.assign(s);
                printf(__FUNCTION__" : de_func %s\n", de_func_name.c_str());
              }
            }
          }
        }

        double val;

        if (!strcmp(s, "position")) { v=fscanf(f, " %lf %lf %lf", &pos()[0], &pos()[1], &pos()[2]); continue; }
        if (!strcmp(s, "direction")) { v=fscanf(f, " %lf %lf %lf", &ahead()[0], &ahead()[1], &ahead()[2]); continue; }
        if (!strcmp(s, "upDirection")) { v=fscanf(f, " %lf %lf %lf", &up()[0], &up()[1], &up()[2]); continue; }

        #define PROCESS(type, name, nameString, doSpline) \
         if (!strcmp(s, nameString)) { v=fscanf(f, " %lf", &val); name = val; continue; }
        PROCESS_CONFIG_PARAMS
        #undef PROCESS

        for (i=0; i<lengthof(par); i++) {
         char p[256];
         sprintf(p, "par%lu", (unsigned long)i);
         if (!strcmp(s, p)) {
           v=fscanf(f, " %f %f %f", &par[i][0], &par[i][1], &par[i][2]);
           break;
         }
        }
      }
      fclose(f);
      printf(__FUNCTION__ " : read '%s'\n", configFile);
      result = true;
     } else {
        printf(__FUNCTION__ " : failed to open '%s'\n", configFile);
     }
     if (result) sanitizeParameters();
     return result;
   }

   // Make sure parameters are OK.
   void sanitizeParameters(void) {
     // Resolution: if only one coordinate is set, keep 4:3 aspect ratio.
     if (width < 1) {
      if (height < 1) { height = 480; }
      width = height*4/3;
     }
     if (height < 1) height = width*3/4;

     // FOV: keep pixels square unless stated otherwise.
     // Default FOV_y is 75 degrees.
     if (fov_x <= 0) {
       if (fov_y <= 0) { fov_y = 75; }
       fov_x = atan(tan(fov_y*PI/180/2)*width/height)/PI*180*2;
     }
     if (fov_y <= 0) fov_y = atan(tan(fov_x*PI/180/2)*height/width)/PI*180*2;

     // Fullscreen: 0=off, anything else=on.
     if (fullscreen != 0 && fullscreen != 1) fullscreen = 1;

     // The others are easy.
     if (multisamples < 1) multisamples = 1;
     //if (speed <= 0) speed = 0.005;  // units/frame
     if (keyb_rot_speed <= 0) keyb_rot_speed = 5;  // degrees/frame
     if (mouse_rot_speed <= 0) mouse_rot_speed = 1;  // degrees/pixel
     if (max_steps < 1) max_steps = 128;
     if (min_dist <= 0) min_dist = 0.0001;
     if (iters < 1) iters = 13;
     if (color_iters < 0) color_iters = 9;
     if (ao_eps <= 0) ao_eps = 0.0005;
     if (ao_strength <= 0) ao_strength = 0.1;
     if (glow_strength <= 0) glow_strength = 0.25;
     if (dist_to_color <= 0) dist_to_color = 0.2;
     if (z_near <= 0) z_near = zNear;
     if (z_far <= 0) z_far = zFar;

     orthogonalize();

     // Don't do anything with user parameters - they must be
     // sanitized (clamped, ...) in the shader.
   }

   // Save configuration.
   void saveConfig(char const* configFile, string* defines = NULL) {
     FILE* f;
     if ((f = fopen(configFile, "w")) != 0) {
       if (defines != NULL)
         fprintf(f, "%s", defines->c_str());
       #define PROCESS(type, name, nameString, doSpline) \
         fprintf(f, nameString " %g\n", (double)name);
       PROCESS_CONFIG_PARAMS
       #undef PROCESS

       fprintf(f, "position %12.12e %12.12e %12.12e\n", pos()[0], pos()[1], pos()[2]);
       fprintf(f, "direction %g %g %g\n", ahead()[0], ahead()[1], ahead()[2]);
       fprintf(f, "upDirection %g %g %g\n", up()[0], up()[1], up()[2]);
       for (size_t i=0; i<lengthof(par); i++) {
         fprintf(f, "par%lu %g %g %g\n", (unsigned long)i, par[i][0], par[i][1], par[i][2]);
       }
       fclose(f);
       printf(__FUNCTION__ " : wrote '%s'\n", configFile);
     }
   }

   // Send parameters to gpu.
   void setUniforms(float x_scale, float x_offset,
                    float y_scale, float y_offset,
                    double speed = 0.0) {
     #define glSetUniformf(name) \
       glUniform1f(glGetUniformLocation(program, #name), name);
     #define glSetUniformfv(name) \
       glUniform3fv(glGetUniformLocation(program, #name), lengthof(name), (float*)name);
     #define glSetUniformi(name) \
       glUniform1i(glGetUniformLocation(program, #name), name);

     glSetUniformf(fov_x); glSetUniformf(fov_y);
     glSetUniformi(max_steps); glSetUniformf(min_dist);
     glSetUniformi(iters); glSetUniformi(color_iters);
     glSetUniformf(ao_eps); glSetUniformf(ao_strength);
     glSetUniformf(glow_strength); glSetUniformf(dist_to_color);
     glSetUniformf(x_scale); glSetUniformf(x_offset);
     glSetUniformf(y_scale); glSetUniformf(y_offset);
     glSetUniformf(speed); glSetUniformi(nrays);
     glSetUniformf(time); glSetUniformf(focus);

     glUniform1f(glGetUniformLocation(program, "xres"), width);

#if defined(PFNGLUNIFORM1DPROC)
     // Pass in double precision values, if supported.
     glUniform1d(glGetUniformLocation(program, "dspeed"), speed);
     glUniform3dv(glGetUniformLocation(program, "deye"), 3, pos());
#endif

     glSetUniformfv(par);

     #undef glSetUniformf
     #undef glSetUniformv
     #undef glUnifor1i
   }

   void render(enum StereoMode stereo) {
     switch(stereo) {
       case ST_OVERUNDER: {  // left / right
         activate();
         setUniforms(1.0, 0.0, 2.0, 1.0, +speed);
         glRects(-1,-1,1,0);  // draw bottom half of screen
         activate();
         setUniforms(1.0, 0.0, 2.0, -1.0, -speed);
         glRects(-1,0,1,1);  // draw top half of screen
         } break;
       case ST_XEYED: {  // right | left
         activate();
         setUniforms(2.0, +1.0, 1.0, 0.0, +speed);
         glRectf(-1,-1,0,1);  // draw left half of screen
         activate();
         setUniforms(2.0, -1.0, 1.0, 0.0, -speed);
         glRectf(0,-1,1,1);  // draw right half of screen
         } break;
       case ST_SIDEBYSIDE: {  // left | right
         activate();
         setUniforms(2.0, +1.0, 1.0, 0.0, -speed);
         glRectf(-1,-1,0,1);  // draw left half of screen
         activate();
         setUniforms(2.0, -1.0, 1.0, 0.0, +speed);
         glRectf(0,-1,1,1);  // draw right half of screen
         } break;
       case ST_NONE:
         activate();
         setUniforms(1.0, 0.0, 1.0, 0.0, speed);
         glRects(-1,-1,1,1);  // draw entire screen
         break;
       case ST_INTERLACED:
         activate();
         setUniforms(1.0, 0.0, 1.0, 0.0, speed);
         glRects(-1,-1,1,1);  // draw entire screen
         break;
      }
    }
} camera,  // Active camera view
  config;  // Global configuration set

vector<KeyFrame> keyframes;  // Keyframes

void suggestDeltaTime(KeyFrame& camera, const vector<KeyFrame>& keyframes) {
  if (keyframes.empty()) {
    camera.delta_time = 0;
  } else {
    double dist = camera.distanceTo(keyframes[keyframes.size() - 1]);
    double steps = dist / camera.speed;
    camera.delta_time = steps / config.fps;
  }
}

#define NSUBFRAMES 1000  // # splined subframes between keyframes.

void CatmullRom(const vector<KeyFrame>& keyframes,
                vector<KeyFrame>* output,
                bool loop = false,
                int nsubframes = NSUBFRAMES) {
  output->clear();
  if (keyframes.size() < 2) return;  // Need at least two points.

  vector<KeyFrame> controlpoints(keyframes);
  if (loop) {
    // Replicate first two at end.
    controlpoints.push_back(keyframes[0]);
    if (controlpoints[controlpoints.size()-1].delta_time == 0) {
      suggestDeltaTime(controlpoints[controlpoints.size() - 1], keyframes);
    }
    controlpoints.push_back(keyframes[1]);
    // Last one is p0 for spline of first keyframe.
    controlpoints.push_back(keyframes[keyframes.size()-1]);
  } else {
    // Replicate last one twice more.
    controlpoints.push_back(keyframes[keyframes.size()-1]);
    controlpoints.push_back(keyframes[keyframes.size()-1]);
    // Last one is p0 for spline of first keyframe.
    controlpoints.push_back(keyframes[0]);
  }

  size_t n = controlpoints.size();

  // Compute controlpoint time based on sum of delta_time up to it.
  // Note we don't spline delta_time but we do spline time.
  float time = 0;
  for (size_t i = 0; i < n - 1; ++i) {
    if (i) time += controlpoints[i].delta_time;
    controlpoints[i].time = time;
  }
  // Last one's time is splined into first one. Set to 0.
  controlpoints[n - 1].time = 0;

  for (size_t i = 0; i < n - 3; ++i) {
    const KeyFrame *p0 = &controlpoints[i>0?i-1:n-1];
    const KeyFrame *p1 = &controlpoints[i];
    const KeyFrame *p2 = &controlpoints[i+1];
    const KeyFrame *p3 = &controlpoints[i+2];
    for (int f = 0; f < nsubframes; ++f) {
      KeyFrame tmp = config;  // Start with default values.
      tmp.setKey(f == 0);
      const double t = f * (1. / nsubframes);

      // The CatmullRom spline function; 0 <= t <= 1
      #define SPLINE(X,p0,p1,p2,p3) \
        (X = (double)(.5 * ( (2 * p1 + \
                            (-p0 + p2)*t + \
                            (2*p0 - 5*p1 + 4*p2 - p3)*t*t + \
                            (-p0 + 3*p1 - 3*p2 + p3)*t*t*t))))

      // Spline position, direction.
      for (size_t j = 0; j < lengthof(tmp.v); ++j) {
        SPLINE(tmp.v[j], p0->v[j], p1->v[j], p2->v[j], p3->v[j]);
      }
      // Spline par[] array.
      for (size_t j = 0; j < lengthof(tmp.par); ++j) {
        SPLINE(tmp.par[j][0],
               p0->par[j][0], p1->par[j][0], p2->par[j][0], p3->par[j][0]);
        SPLINE(tmp.par[j][1],
               p0->par[j][1], p1->par[j][1], p2->par[j][1], p3->par[j][1]);
        SPLINE(tmp.par[j][2],
               p0->par[j][2], p1->par[j][2], p2->par[j][2], p3->par[j][2]);
      }
      // Spline all other params. Some are non-sensical.
      #define PROCESS(a,b,c,doSpline) \
        if (doSpline) SPLINE(tmp.b, p0->b, p1->b, p2->b, p3->b);
      PROCESS_CONFIG_PARAMS
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
  CTL_FOV = lengthof(camera.par), CTL_RAY, CTL_ITER, CTL_AO, CTL_GLOW,
  CTL_TIME, CTL_CAM, CTL_3D,
  CTL_LAST = CTL_3D,
} Controller;


// Controller modifiers.
// They get a pointer to the modified value and a signed count of consecutive changes.
void m_mul(float* x, int d) { *x *= pow(10, sign(d)/20.); }
void m_mulSlow(float* x, int d) { *x *= pow(10, sign(d)/40.); }
void m_mulSlow(double* x, int d) { *x *= pow(10, sign(d)/40.); }
void m_tan(float* x, int d) { *x = atan(tan(*x*PI/180/2) * pow(0.1, sign(d)/40.) ) /PI*180*2; }
void m_progressiveInc(int* x, int d) { *x += sign(d) * ((abs(d)+4) / 4); }
void m_progressiveAdd(float* x, int d) { *x += 0.001 * (sign(d) * ((abs(d)+4) / 4)); }
void m_singlePress(int* x, int d) { if (d==1 || d==-1) *x += d; }

void m_rotateX(int d) { camera.rotate(sign(d)*camera.keyb_rot_speed, 0, 1, 0); }
void m_rotateY(int d) { camera.rotate(-sign(d)*camera.keyb_rot_speed, 1, 0, 0); }

void m_rotateX2(float d) { camera.rotate(d, 0, 1, 0); }
void m_rotateY2(float d) { camera.rotate(d, 1, 0, 0);}
void m_rotateZ2(float d) { camera.rotate(d, 0, 0, 1); }

// Print controller values into a string.
char* printController(char* s, Controller c) {
  assert(c <= CTL_LAST);
  switch (c) {
    default: {
      char x[8],y[8]; sprintf(x, "par%dx", c); sprintf(y, "par%dy", c);
      sprintf(s, "%s %.3f %s %.3f",
        parName[c][1] ? parName[c][1] : y, camera.par[c][1],
        parName[c][0] ? parName[c][0] : x, camera.par[c][0]
      );
    } break;
    case CTL_FOV: sprintf(s, "Fov %.3g %.3g", camera.fov_x, camera.fov_y); break;
    case CTL_RAY: sprintf(s, "Ray %.2e steps %d", camera.min_dist, camera.max_steps); break;
    case CTL_ITER: sprintf(s, "It %d|%d", camera.iters, camera.color_iters); break;
    case CTL_AO: sprintf(s, "aO %.2e aOeps %.2e", camera.ao_strength, camera.ao_eps); break;
    case CTL_GLOW: sprintf(s, "Glow %.3f Dist %.2e", camera.glow_strength, camera.dist_to_color); break;
    case CTL_CAM: {
      sprintf(s, "Look [%4d %4d %4d]",
              (int)(camera.ahead()[0]*100),
              (int)(camera.ahead()[1]*100),
              (int)(camera.ahead()[2]*100));
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
    default: m_progressiveAdd(&camera.par[c][1], d); break;
    case CTL_FOV: m_tan(&camera.fov_x, d); m_tan(&camera.fov_y, d); break;
    case CTL_RAY: m_mul(&camera.min_dist, d); break;
    case CTL_ITER: m_singlePress(&camera.iters, d); if ((camera.iters&1)==(d>0)) m_singlePress(&camera.color_iters, d); break;
    case CTL_AO: m_mulSlow(&camera.ao_strength, d); break;
    case CTL_GLOW: m_progressiveAdd(&camera.glow_strength, d); break;
    case CTL_CAM: m_rotateY(d); break;
    case CTL_TIME: m_progressiveAdd(&camera.delta_time, d); break;
    case CTL_3D: m_progressiveAdd(&camera.focus, d * 100); break;
  }
  // Enforce sane bounds.
  if (camera.delta_time < 0) camera.delta_time = 0;
}

// Update controller.x by the signed count of consecutive changes.
void updateControllerX(Controller c, int d, bool alt) {
  assert(c <= CTL_LAST);
  switch (c) {
    default: m_progressiveAdd(&camera.par[c][alt?2:0], d); break;
    case CTL_FOV: m_tan(&camera.fov_x, d); break;
    case CTL_RAY: m_progressiveInc(&camera.max_steps, d); break;
    case CTL_ITER: m_singlePress(&camera.color_iters, d); break;
    case CTL_AO: m_mul(&camera.ao_eps, d); break;
    case CTL_GLOW: m_mul(&camera.dist_to_color, d); break;
    case CTL_CAM: m_rotateX(d); break;
    case CTL_TIME: m_mulSlow(&camera.speed, d); break;
    case CTL_3D: m_mulSlow(&camera.speed, d); break;
  }
}

// Change the active controller with a keypress.
void changeController(SDLKey key, Controller* c) {
  switch (key) {
    case SDLK_0: case SDLK_KP0: *c = (Controller)0; break;
    case SDLK_1: case SDLK_KP1: *c = (Controller)1; break;
    case SDLK_2: case SDLK_KP2: *c = (Controller)2; break;
    case SDLK_3: case SDLK_KP3: *c = (Controller)3; break;
    case SDLK_4: case SDLK_KP4: *c = (Controller)4; break;
    case SDLK_5: case SDLK_KP5: *c = (Controller)5; break;
    case SDLK_6: case SDLK_KP6: *c = (Controller)6; break;
    case SDLK_7: case SDLK_KP7: *c = (Controller)7; break;
    case SDLK_8: case SDLK_KP8: *c = (Controller)8; break;
    case SDLK_9: case SDLK_KP9: *c = (Controller)9; break;
    case SDLK_f: *c = CTL_FOV; break;
    case SDLK_g: *c = CTL_GLOW; break;
    case SDLK_l: *c = CTL_CAM; break;
    case SDLK_i: *c = CTL_ITER; break;
    case SDLK_o: *c = CTL_AO; break;
    case SDLK_r: *c = CTL_RAY; break;
    case SDLK_t: *c = CTL_TIME; break;
    case SDLK_v: *c = CTL_3D; break;
    default: break;  // no change
  }
}


////////////////////////////////////////////////////////////////
// Graphics.

// Position of the OpenGL window on the screen.
int viewportOffset[2];

// Is the mouse and keyboard input grabbed?
int grabbedInput = 1;

void saveScreenshot(char const* tgaFile) {
  TGA tga;
  tga.readFramebuffer(config.width, config.height, viewportOffset);
  if (tga.writeFile(tgaFile))
    printf(__FUNCTION__ " : wrote %s\n", tgaFile);
  else
    printf(__FUNCTION__ " : failed to write %s\n", tgaFile);
}

TGA background;

void LoadBackground() {
  background.readFile("background.tga");
  if (background.data()) {
    printf(__FUNCTION__ " : loaded background image from '%s'\n", "background.tga");
  }
}

// Return BGR value of pixel x,y.
unsigned int getBGRpixel(int x, int y) {
  unsigned char img[3];
  int height = config.height;
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_FRONT);
  glReadPixels(viewportOffset[0] + x, viewportOffset[1] + height - 1 - y, 1, 1, GL_BGR, GL_UNSIGNED_BYTE, img);
  unsigned int val =
      img[0] * 256 * 256 +
      img[1] * 256 +
      img[2];
  return val;
}

string glsl_source;

// Compile and activate shader programs. Return the program handle.
int setupShaders(void) {
  char const* vs;
  char const* fs;
  GLuint v,f,p;
  char log[2048]; int logLength;

  (vs = readFile(VERTEX_SHADER_FILE)) || ( vs = default_vs );
  (fs = readFile(FRAGMENT_SHADER_FILE)) || ( fs = default_fs );

  if (fs != default_fs) {
     printf(__FUNCTION__ " : read shader from %s\n", FRAGMENT_SHADER_FILE);
  } else {
     printf(__FUNCTION__ " : using default shader\n");
  }
  if (vs != default_vs) {
     printf(__FUNCTION__ " : read vertex shader from %s\n", VERTEX_SHADER_FILE);
  } else {
     printf(__FUNCTION__ " : using default vertex shader\n");
  }

  glsl_source.assign(defines + fs);

  p = glCreateProgram();

  v = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(v, 1, &vs, 0);
  glCompileShader(v);
  glGetShaderInfoLog(v, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  f = glCreateShader(GL_FRAGMENT_SHADER);
  if (!defines.empty()) {
    const char* srcs[2] = {defines.c_str(), fs};
    glShaderSource(f, 2, srcs, 0);
  } else {
    glShaderSource(f, 1, &fs, 0);
  }
  glCompileShader(f);
  glGetShaderInfoLog(f, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  glAttachShader(p, v);
  glAttachShader(p, f);
  glLinkProgram(p);

  glGetProgramInfoLog(p, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  if (vs != default_vs) free((char*)vs);
  if (fs != default_fs) free((char*)fs);

  GLint status;
  glGetProgramiv(p, GL_LINK_STATUS, &status);
  if (status != GL_TRUE) die("setupShaders() fails\n");

  return p>0?p:0;
}

// Compile and activate shader programs for frame buffer manipulation.
// Return the program handle.
int setupShaders2(void) {
  char const* vs;
  char const* fs;
  GLuint v,f,p;
  char log[2048]; int logLength;

  (vs = readFile(FRAME_VERTEX_SHADER_FILE)) || ( vs = frame_default_vs );
  (fs = readFile(FRAME_FRAGMENT_SHADER_FILE)) || ( fs = frame_default_fs );

  if (fs != frame_default_fs) {
     printf(__FUNCTION__ " : read shader from %s\n", FRAME_FRAGMENT_SHADER_FILE);
  } else {
     printf(__FUNCTION__ " : using default shader\n");
  }

  glsl_source.append(fs);

  p = glCreateProgram();

  v = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(v, 1, &vs, 0);
  glCompileShader(v);
  glGetShaderInfoLog(v, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  f = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(f, 1, &fs, 0);
  glCompileShader(f);
  glGetShaderInfoLog(f, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  glAttachShader(p, v);
  glAttachShader(p, f);
  glLinkProgram(p);

  glGetProgramInfoLog(p, sizeof(log), &logLength, log);
  if (logLength) fprintf(stderr, __FUNCTION__ " : %s\n", log);

  if (vs != frame_default_vs) free((char*)vs);
  if (fs != frame_default_fs) free((char*)fs);

  GLint status;
  glGetProgramiv(p, GL_LINK_STATUS, &status);
  if (status != GL_TRUE) die("setupShaders2() fails\n");

  return p>0?p:0;
}

// Detach & delete any shaders, delete program.
void cleanupShaders(int p) {
  GLuint shaders[2];
  GLsizei count = 2;
  if (!glIsProgram(p)) return;
  glGetAttachedShaders(p, count, &count, shaders);
  for (GLsizei i = 0; i < count; ++i) {
    glDetachShader(p, shaders[i]);
    glDeleteShader(shaders[i]);
  }
  glDeleteProgram(p);
}

void changeWorkingDirectory(const char* configFile) {
  char dirName[MAX_PATH];
  strncpy(dirName, configFile, sizeof dirName);
  dirName[sizeof dirName - 1] = 0;
  if (strstr(configFile, ".cfg.data") == NULL) {
    strncat(dirName, ".data", sizeof dirName);
    MKDIR(dirName);
  } else {
    // Already have a .cfg.data in path, chdir to it.
    int i = strlen(dirName);
    while (i > 0 && strchr("/\\", dirName[--i]) == NULL);
    dirName[i] = 0;
  }
  if (!chdir(dirName)) printf(__FUNCTION__ " : chdir '%s'\n", dirName);
}

// Initializes the video mode, OpenGL state, shaders, camera and shader parameters.
// Exits the program if an error occurs.
void initGraphics() {
  GLenum status = GL_NO_ERROR;
  // If not fullscreen, use the color depth of the current video mode.
  int bpp = 24;  // FSAA works reliably only in 24bit modes
  if (!config.fullscreen) {
    SDL_VideoInfo const* info = SDL_GetVideoInfo();
    bpp = info->vfmt->BitsPerPixel;
  }

  // Set attributes for the OpenGL window.
  SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, config.depth_size);

  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  if (config.multisamples == 1) {
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 0);
  } else {
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, config.multisamples);
  }

  // Set the video mode, hide the mouse and grab keyboard and mouse input.
  if (screen == NULL) SDL_putenv((char*)"SDL_VIDEO_CENTERED=center");
  if (screen != NULL) SDL_FreeSurface(screen);
  (screen = SDL_SetVideoMode(config.width, config.height, bpp,
      SDL_OPENGL | (config.fullscreen ? SDL_FULLSCREEN : SDL_RESIZABLE)))
    || die("Video mode initialization failed: %s\n", SDL_GetError());

  if (config.multisamples > 1) {
     glEnable(GL_MULTISAMPLE);  // redundant?
     int i;
     SDL_GL_GetAttribute(SDL_GL_MULTISAMPLESAMPLES, &i);
     printf(__FUNCTION__ " : multi says %d, asked for %d\n",
            i, config.multisamples);
  }

  SDL_GL_GetAttribute(SDL_GL_DEPTH_SIZE, &config.depth_size);
  printf(__FUNCTION__ " : depth size %u\n", config.depth_size);

  if (grabbedInput) {
    // order is important
    SDL_ShowCursor(SDL_DISABLE);
    SDL_WM_GrabInput(SDL_GRAB_ON);
  } else {
    SDL_ShowCursor(SDL_ENABLE);
    SDL_WM_GrabInput(SDL_GRAB_OFF);
  }

  // TODO: logging.
  //int samples = 0; glGetIntegerv(GL_SAMPLES, &samples); printf("%dx%d, %d bpp, FSAA %d\n", screen->w, screen->h, screen->format->BitsPerPixel, samples);

  // If we got higher resolution (which can happen in fullscreen mode),
  // use a centered viewport.
  viewportOffset[0] = (screen->w - config.width)/2;
  viewportOffset[1] = (screen->h - config.height)/2;
  glViewport(viewportOffset[0], viewportOffset[1], config.width, config.height);

  if (hand_cursor == NULL) hand_cursor = init_system_cursor(kHand);
  if (arrow_cursor == NULL) arrow_cursor = SDL_GetCursor();

  // Enable shader functions and compile shaders.
  // Needs to be done after setting the video mode.
  enableShaderProcs() ||
      die("This program needs support for GLSL shaders.\n");

  cleanupShaders(program);
  cleanupShaders(dof_program);

  (program = setupShaders()) ||
      die("Error in GLSL shader compilation (see stderr.txt for details).\n");

  if (background_texture == 0 && background.data() != NULL) {
    // Load background image into texture
    glGenTextures(1, &background_texture);
    glBindTexture(GL_TEXTURE_2D, background_texture);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, background.width(), background.height(),
                 0, GL_BGR, GL_UNSIGNED_BYTE, background.data());
    glGenerateMipmap(GL_TEXTURE_2D);
    printf(__FUNCTION__ " : background texture at %d\n", background_texture);
    glBindTexture(GL_TEXTURE_2D, 0);
  }

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

  if (!config.enable_dof) return;

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

  (dof_program = setupShaders2()) ||
      die("Error in GLSL shader compilation (see stderr.txt for details).\n");


  // Create depthbuffer
  glDeleteRenderbuffers(1, &depthBuffer);
  glGenRenderbuffers(1, &depthBuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                        config.width, config.height);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

  // Create texture to render to
  glDeleteTextures(1, &texture);
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);

  // Allocate storage
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, config.width, config.height,
               0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);

  glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

  // Allocate / generate mips
  glGenerateMipmap(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, 0);

  // Create framebuffer
  glDeleteFramebuffers(1, &fbo);
  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  // Attach texture to framebuffer as colorbuffer
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                         GL_TEXTURE_2D, texture, 0);

  // Attach depthbuffer to framebuffer
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, depthBuffer);

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

  status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
  if (status != GL_FRAMEBUFFER_COMPLETE)
    die(__FUNCTION__ " : glCheckFramebufferStatus() : %04x\n", status);

  // Back to normal framebuffer.
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);
}

TwBar* bar = NULL;

// Find '\n#define foo par[x].z  // {twbar params}' in glsl_source.
void initTwParDefines() {
  size_t start = 0;
  while ((start = glsl_source.find("\n#define ", start + 1)) != string::npos) {
    size_t eol = glsl_source.find("\n", start + 1);
    if (eol == string::npos) continue;
    string line(glsl_source, start + 1, eol - (start + 1));

    size_t parStart = line.find(" par[");
    if (parStart == string::npos || parStart < 8) continue;
    string varName(line, 8, parStart - 8);

    size_t attr_start = line.find("{");
    size_t attr_end = line.find("}");
    string attr;
    if (attr_start != string::npos &&
        attr_end != string::npos &&
        attr_end > attr_start)
      attr.assign(line, attr_start + 1, attr_end - (attr_start + 1));

    int index;
    char _xyz = 'x';
    if (sscanf(line.c_str() + parStart + 5, "%d].%c",
               &index, &_xyz) < 1) continue;
    if (index < 0 || index > (int)lengthof(camera.par)) continue;
    if (_xyz < 'x' || _xyz > 'z') continue;

    printf("parameter %s par[%d].%c {%s}\n",
           varName.c_str(), index, _xyz, attr.c_str());

    float* address = &camera.par[index][_xyz-'x'];

    if (varName.find("Color") != string::npos) {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_COLOR3F, address, "");
    } else if (varName.find("Vector") != string::npos) {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_DIR3F, address, "");
    } else {
      TwAddVarRW(bar, varName.c_str(), TW_TYPE_FLOAT, address, attr.c_str());
    }
  }
}

// Find '\nuniform <type> name;' in glsl_source.
// Pick up type and twbar params, if found.
void initTwUniform(const string& name, void* addr) {
  size_t start = 0;
  while ((start = glsl_source.find("\nuniform ", start + 1)) != string::npos) {
    size_t eol = glsl_source.find("\n", start + 1);
    string line(glsl_source, start + 1, eol - (start + 1));

    size_t attr_start = line.find("{");
    size_t attr_end = line.find("}");
    string attr;
    if (attr_start != string::npos &&
        attr_end != string::npos &&
        attr_end > attr_start)
      attr.assign(line, attr_start + 1, attr_end - (attr_start + 1));

    if (line.find("uniform float " + name + ";") == 0) {
      TwAddVarRW(bar, name.c_str(), TW_TYPE_FLOAT, (float*)addr, attr.c_str());
    } else if (line.find("uniform int " + name + ";") == 0) {
      TwAddVarRW(bar, name.c_str(), TW_TYPE_INT32, (int*)addr, attr.c_str());
    }
  }
}

void initTwBar() {
  if (bar == NULL) TwInit(TW_OPENGL, NULL);

  TwWindowSize(config.width, config.height);

  if (bar != NULL) return;

  bar = TwNewBar("boxplorer");

  #define PROCESS(a,b,c,d) initTwUniform(c, &camera.b);
  PROCESS_CONFIG_PARAMS
  #undef PROCESS

  initTwParDefines();
}

void LoadKeyFrames(bool fixedFov) {
  char filename[256];
  for (int i = 0; ; ++i) {
    sprintf(filename, "%s-%u.cfg", kKEYFRAME, i);
    KeyFrame tmp;
    if (!tmp.loadConfig(filename)) break;
    if (fixedFov) {
      tmp.width = config.width;
      tmp.height = config.height;
      tmp.fov_x = config.fov_x;
      tmp.fov_y = config.fov_y;
    }
    keyframes.push_back(tmp);
  }
  printf(__FUNCTION__ " : loaded %lu keyframes\n",
      (unsigned long)keyframes.size());
}

void SaveKeyFrames() {
  char filename[256];
  for (size_t i = 0; i < keyframes.size(); ++i) {
    sprintf(filename, "%s-%lu.cfg", kKEYFRAME, (unsigned long)i);
    keyframes[i].saveConfig(filename);
  }
}

////////////////////////////////////////////////////////////////
// Setup, input handling and drawing.

int main(int argc, char **argv) {
  bool rendering = false;
  bool loop = false;
  bool useTime = false;
  bool configSpeed = false;
  bool fixedFov = false;
  int enableDof = 0;
#if defined(_MACOSX)
  // Shit won't work on mac.
  enableDof = -1;
#endif

  // Peel known options off the back..
  while (argc>1) {
    if (!strcmp(argv[argc-1], "--overunder")) {
      stereoMode = ST_OVERUNDER;
    } else if (!strcmp(argv[argc-1], "--interlaced")) {
      stereoMode = ST_INTERLACED;
      defines.append("#define ST_INTERLACED\n");
    } else if (!strcmp(argv[argc-1], "--xeyed")) {
      stereoMode = ST_XEYED;
    } else if (!strcmp(argv[argc-1], "--sidebyside")) {
      stereoMode = ST_SIDEBYSIDE;
    } else if (!strcmp(argv[argc-1], "--render")) {
      rendering = true;
    } else if (!strcmp(argv[argc-1], "--time")) {
      useTime = true;
    } else if (!strcmp(argv[argc-1], "--speed")) {
      configSpeed = true;
    } else if (!strcmp(argv[argc-1], "--disabledof")) {
      enableDof = -1;
    } else if (!strcmp(argv[argc-1], "--enabledof")) {
      enableDof = 1;
    } else if (!strcmp(argv[argc-1], "--fixedfov")) {
      fixedFov = true;
    } else if (!strcmp(argv[argc-1], "--loop")) {
      loop = true;
    } else if (!strncmp(argv[argc-1], "--kf=", 5)) {
      kKEYFRAME = argv[argc-1] + 5;
    } else break;
    --argc;
  }

  if (stereoMode == ST_NONE) defines.append("#define ST_NONE\n");

  const char* configFile = (argc>=2 ? argv[1] : DEFAULT_CONFIG_FILE);

  // Load configuration.
  if (config.loadConfig(configFile, &defines)) {
    changeWorkingDirectory(configFile);
  } else {
    die("Usage: boxplorer <configuration-file.cfg>\n");
  }

  // Sanitize / override config parameters.
  if (loop) config.loop = true;
  if (enableDof) config.enable_dof = (enableDof == 1);  // override
  if (stereoMode == ST_INTERLACED) config.enable_dof = 0;  // mipmapping does not work for interlaced.
  if (config.fps < 5) config.fps = 30;
  if (config.depth_size < 16) config.depth_size = 16;
  if (stereoMode == ST_XEYED) config.width *= 2;

  camera = config;

  bool keyframesChanged = false;
  LoadKeyFrames(fixedFov);
  LoadBackground();
 
  // Initialize SDL and OpenGL graphics.
  SDL_Init(SDL_INIT_VIDEO) == 0 ||
      die("SDL initialization failed: %s\n", SDL_GetError());
  atexit(SDL_Quit);
  
   // Set up the video mode, OpenGL state, shaders and shader parameters.
  initGraphics();
  initTwBar();
  initFPS(FPS_FRAMES_TO_AVERAGE);

  printf(__FUNCTION__ " : GL_EXTENSIONS: %s\n", glGetString(GL_EXTENSIONS));

  // Main loop.
  Controller ctl = CTL_CAM;  // the default controller is camera rotation
  int consecutiveChanges = 0;

  int done = 0;
  int frameno = 0;

  vector<KeyFrame> splines;
  size_t splines_index = 0;

  float frame_time = 1 / config.fps;
  float render_time = 0;
  float render_start = 0;

  float dist_along_spline = 0;
  size_t keyframe = keyframes.size();

  bool ignoreNextMouseUp = false;
  bool dragging = false;

  if (rendering) {
    // Rendering a sequence to disk. Spline the keyframes.
    CatmullRom(keyframes, &splines, config.loop);
  }

  double last_de = 0;

  // Check availability of DE; setup func ptr.
  if (!de_func_name.empty()) {
    if (DE64_funcs.find(de_func_name) != DE64_funcs.end()) {
      // First pick is double version.
      de_func_64 = DE64_funcs[de_func_name];
    } else if (DE_funcs.find(de_func_name) != DE_funcs.end()) {
      // Next pick float version.
      de_func = DE_funcs[de_func_name];
    } else {
      printf(__FUNCTION__ ": unknown DE %s\n", de_func_name.c_str());
      de_func_name.clear();
    }
  }

  while (!done) {
    int ctlXChanged = 0, ctlYChanged = 0;

    // Splined keyframes playback logic. Mess.
    if (!splines.empty()) {
      if (rendering && splines_index >= splines.size()) break;  // done

      // Loop if asked to.
      if (config.loop && splines_index >= splines.size()) {
        splines_index = 0;
        render_time = 0;
        render_start = now();
      }

      // Figure out whether to draw a splined frame or skip to next.
      if (splines_index < splines.size()) {
        if (!config.loop && splines_index == splines.size() - 1) {
          camera = keyframes[keyframes.size() - 1];  // End at last position.
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
            if (camera.time < render_time) continue;
            render_time += frame_time;
          } else {
            // Previewing. Use real time (low framerate == jumpy preview!).
            float n = now();
            if (n > render_start + camera.time) continue;  // late, skip frame.
            float w = (render_start + camera.time) - n;
            if (w >= frame_time) {  // early, redraw frame.
              splines_index = prev_splines_index;
            }
          }
        } else {
         if (dist_along_spline < (configSpeed?config.speed:camera.speed))
           continue;
        }
        dist_along_spline -= (configSpeed?config.speed:camera.speed);
      } else {
        splines.clear();
      }
    }

    if (config.enable_dof && camera.dof_scale > .0001) {
      // If we have some DoF to render, render to texture.
      glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    }

    if (de_func || de_func_64) {
      // Get a DE.
      // Setup just the vars needed for DE. For now, iters and par[0..10]
      // TODO: make de() method of camera?
      GLSL::iters = camera.iters;
      for (int i = 0; i < 10; ++i) {
        GLSL::par[i] = GLSL::vec3(camera.par[i][0], camera.par[i][1], camera.par[i][2]);
      }
      double de =
        de_func_64
          ?GLSL::abs(de_func_64(
              GLSL::dvec3(camera.pos()[0], camera.pos()[1], camera.pos()[2])))
          :GLSL::abs(de_func(
              GLSL::vec3(camera.pos()[0], camera.pos()[1], camera.pos()[2])));
      if (de != last_de) {
        printf("de=%12.12e\n", de);
        camera.speed = GLSL::clamp(de/10.0,DBL_EPSILON,1.0);
      }
      last_de = de;
    }

    if (!rendering) {
      // If we're rendering a sequence to disk, we don't care about z-buffer.
      // Otherwise, just overwrite since we write every pixel.
      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_ALWAYS);
    }

    if (background_texture) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, background_texture);
      glUseProgram(program);
      glUniform1i(glGetUniformLocation(program, "bg_texture"), 0);
    } else {
      glUseProgram(program);
    }
    glUniform1i(glGetUniformLocation(program, "use_bg_texture"), background_texture);

    camera.render(stereoMode);

    glUseProgram(0);

    if (background_texture) {
      glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!config.no_spline &&
        stereoMode == ST_NONE &&
        keyframes.size() > 1 &&
        splines.empty()) {
      // Draw keyframe splined path, if we have 2+ keyframes and not rendering
      glDepthFunc(GL_LESS);

      camera.activateGl();

      vector<KeyFrame> splines;
      CatmullRom(keyframes, &splines, config.loop);

      glColor4f(.8,.8,.8,1);
      glLineWidth(1);
      glEnable(GL_LINE_SMOOTH);

      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // right eye
      for (size_t i = 0; i < splines.size(); ++i) {
        splines[i].move(splines[i].speed, 0, 0);
        glVertex3dv(splines[i].pos());
        splines[i].move(-2*splines[i].speed, 0, 0);
      }
      glEnd();
      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // left eye
      for (size_t i = 0; i < splines.size(); ++i) {
        glVertex3dv(splines[i].pos());
      }
      glEnd();

#if 0
      // draw light 1 path
      glColor4f(.1,.9,.9,1);
      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // light
      for (size_t i = 0; i < splines.size(); ++i) {
        glVertex3fv(splines[i].par[19]);
      }
      glEnd();
#endif

      glLineWidth(13);
      for (size_t i = 0; i < keyframes.size(); ++i) {
        glColor4f(0,0,1 - (i/256.0),1);  // Encode keyframe # in color.
        glBegin(GL_LINES);
        KeyFrame tmp = keyframes[i];
        tmp.move(tmp.speed, 0, 0);
        glVertex3dv(tmp.pos());
        tmp = keyframes[i];
        tmp.move(-tmp.speed, 0, 0);
        glVertex3dv(tmp.pos());
        glEnd();
      }
    }

    glDisable(GL_DEPTH_TEST);

    if (config.enable_dof && camera.dof_scale > .0001) {
      // If we're rendering some DoF, draw texture on screen.
      glBindFramebuffer(GL_FRAMEBUFFER, 0);

      glEnable(GL_TEXTURE_2D);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture);
      glGenerateMipmap(GL_TEXTURE_2D);

      glUseProgram(dof_program);  // Activate our alpha channel shader.

      glUniform1i(glGetUniformLocation(dof_program, "my_texture"), 0);
      glUniform1f(glGetUniformLocation(dof_program, "z_near"), camera.z_near);
      glUniform1f(glGetUniformLocation(dof_program, "z_far"), camera.z_far);
      glUniform1f(glGetUniformLocation(dof_program, "dof_scale"),
                  camera.dof_scale);
      glUniform1f(glGetUniformLocation(dof_program, "dof_offset"),
                  camera.dof_offset);

      // Ortho projection, entire screen in regular pixel coordinates.
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(0, config.width, config.height, 0, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      // Draw our texture covering entire screen, running the frame shader.
      glColor4f(1,1,1,1);
      glBegin(GL_QUADS);
        glTexCoord2f(0,1);
        glVertex2f(0,0);
        glTexCoord2f(0,0);
        glVertex2f(0,config.height);
        glTexCoord2f(1,0);
        glVertex2f(config.width,config.height);
        glTexCoord2f(1,1);
        glVertex2f(config.width,0);
      glEnd();

      glUseProgram(0);

      glDisable(GL_TEXTURE_2D);
    }

    // Draw AntTweakBar
    if (!grabbedInput && !rendering) TwDraw();

    {
      GLenum err = glGetError();
      if (err != GL_NO_ERROR) printf("glGetError():%04x\n", err);
    }

    SDL_GL_SwapBuffers();

    if (rendering && !splines.empty()) {
      // If we're playing a sequence back, save every frame to disk.
      char filename[256];
      sprintf(filename, "frame-%05d.tga", frameno);
      ++frameno;
      saveScreenshot(filename);
    }

    updateFPS();

    // Show position and fps in the caption.
    char caption[2048], controllerStr[256];
    sprintf(caption, "%s %.2ffps %5lu [%.3lf %.3lf %.3lf] %dms",
        printController(controllerStr, ctl),
        getFPS(), (unsigned long)splines_index,
        camera.pos()[0], camera.pos()[1], camera.pos()[2],
        getLastFrameDuration()
    );
    SDL_WM_SetCaption(caption, 0);

    // Process events.
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (grabbedInput ||
          !TwEventSDL(&event, SDL_MAJOR_VERSION, SDL_MINOR_VERSION))
      switch (event.type) {
      case SDL_VIDEORESIZE: {
            config.width = event.resize.w;
            config.height = event.resize.h;
            config.fov_x = 0;  // go for square pixels..
            config.sanitizeParameters();
            grabbedInput = 1; initGraphics(); initTwBar();
            printf("resize(%d, %d)\n", config.width, config.height);
      } break;

      case SDL_QUIT: done |= 1; break;

      case SDL_MOUSEBUTTONDOWN: {
        if (grabbedInput == 0) {
          unsigned int bgr = getBGRpixel(event.button.x, event.button.y);
          if ((bgr & 0xffff) == 0) {
            // No red or green at all : probably a keyframe marker (fragile).
            size_t kf = 255 - (bgr >> 16);
            if (kf < keyframes.size()) {
               printf("selected keyframe %lu\n", (unsigned long)kf);
               keyframe = kf;
               dragging = true;
            }
          }
        } else switch(event.button.button) {
          case 1:  // left mouse
            grabbedInput = 0;
            SDL_ShowCursor(SDL_ENABLE);
            SDL_WM_GrabInput(SDL_GRAB_OFF);
            ignoreNextMouseUp = true;
            break;
          case 4:  // mouse wheel up, increase eye distance
             camera.speed *= 1.1;
            break;
          case 5:  // mouse wheel down, decrease eye distance
            camera.speed *= 0.9;
            break;
        }
      } break;

      case SDL_MOUSEMOTION: {
         if (grabbedInput == 0) {
           if (!dragging) {
              unsigned int bgr = getBGRpixel(event.motion.x, event.motion.y);
              if ((bgr & 0xffff) == 0) {
                // No red or green at all : probably a keyframe marker.
                size_t kf = 255 -(bgr >> 16);
                if (kf < keyframes.size()) {
                  printf("keyframe %lu\n", (unsigned long)kf);
                  SDL_SetCursor(hand_cursor);
                } else {
                  SDL_SetCursor(arrow_cursor);
                }
              } else {
                 SDL_SetCursor(arrow_cursor);
              }
           } else {
             // Drag currently selected keyframe around.
             if (keyframe < keyframes.size()) {
               // Should really be some screenspace conversion..
               // but this works ok for now.
               keyframes[keyframe].moveAbsolute(camera.up(),
                   event.motion.yrel*-.05*keyframes[keyframe].speed);
               keyframes[keyframe].moveAbsolute(camera.right(),
                   event.motion.xrel*.05*keyframes[keyframe].speed);
             }
           }
         }
      } break;

      case SDL_MOUSEBUTTONUP: {
         if (ignoreNextMouseUp == false && grabbedInput == 0) {
           grabbedInput = 1;
           SDL_SetCursor(arrow_cursor);
           SDL_ShowCursor(SDL_DISABLE);
           SDL_WM_GrabInput(SDL_GRAB_ON);
         }
         ignoreNextMouseUp = false;
         dragging = false;
      } break;

      case SDL_KEYDOWN: {
      bool hasAlt = event.key.keysym.mod & (KMOD_LALT|KMOD_RALT);
      bool hasCtrl = event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL);

      switch (event.key.keysym.sym) {
      case SDLK_ESCAPE: {
         if (grabbedInput && !config.fullscreen) {
           grabbedInput = 0;
           SDL_ShowCursor(SDL_ENABLE);
           SDL_WM_GrabInput(SDL_GRAB_OFF);
         } else done |= 1;
      } break;

      // Switch fullscreen mode (loses the whole OpenGL context in Windows).
      case SDLK_RETURN: case SDLK_KP_ENTER: {
        config.fullscreen ^= 1; grabbedInput = 1;
        if (config.fullscreen) {
          if (stereoMode == ST_INTERLACED) {
            // Zalman
            config.height = 1080; config.width = 1920;
          } else if (stereoMode == ST_SIDEBYSIDE || stereoMode == ST_OVERUNDER) {
            // HMZ-T1
            config.height = 720; config.width = 1280;
          }
        }
        initGraphics(); initTwBar();
      } break;

      // Save config and screenshot (filename = current time).
      case SDLK_SYSREQ:
      case SDLK_PRINT: {
        time_t t = time(0);
        struct tm* ptm = localtime(&t);
        char filename[256];
        strftime(filename, 256, "%Y%m%d_%H%M%S.cfg", ptm);
        camera.saveConfig(filename);
        strftime(filename, 256, "%Y%m%d_%H%M%S.tga", ptm);
        saveScreenshot(filename);
      } break;

      // Splined path control.
      case SDLK_HOME: {  // Start spline playback from start.
        if (!keyframes.empty()) {
            CatmullRom(keyframes, &splines, config.loop);
            splines_index = 0;
            render_time = 0;
            render_start = now();
            keyframe = 0;
        } else camera = config;
      } break;

      case SDLK_END: {  // Stop spline playback.
        splines.clear();
        if (!keyframes.empty()) {
          keyframe = keyframes.size() - 1;
          camera = keyframes[keyframe];
        } else camera = config;
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
              camera = keyframes[keyframe];
          }
        }
      } break;

      case SDLK_SPACE:
      case SDLK_INSERT: {  // Add keyframe.
        splines.clear();
        size_t index = keyframes.size();
        if (hasCtrl) {
          // Replace currently selected keyframe.
          index = keyframe;
        }
        if (!hasCtrl) {
          // Need an estimate for delta_time for this new frame.
          suggestDeltaTime(camera, keyframes);
          keyframes.push_back(camera);  // Add keyframe at end.
        } else {
          keyframes[index] = camera;  // Overwrite current keyframe.
        }
        char filename[256];
        sprintf(filename, "%s-%lu.cfg", kKEYFRAME, (unsigned long)index);
        camera.saveConfig(filename);
      } break;

      case SDLK_PAGEUP: {  // Jump keyframe back in spline playback.
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
          splines_index = splines.size()-1;
        }
      } break;

      case SDLK_TAB: {  // Step thru keyframes.
        if (!splines.empty()) {
           // find most recently played keyframe
           int nkeys = 0;
           for (size_t i = 0; i < splines_index; ++i)
              if (splines[i].isKey()) ++nkeys;
           keyframe = nkeys;
        }
        if (!hasCtrl) ++keyframe;
        if (keyframe >= keyframes.size()) keyframe = 0;

        if (keyframe < keyframes.size() && splines.empty()) {
           // Don't jump camera ahead if we were playing, just stop in place.
           camera = keyframes[keyframe];
        }

        if (keyframe < keyframes.size()) {
          printf("at keyframe %lu, speed %.8e, delta_time %f\n",
              (unsigned long)keyframe, keyframes[keyframe].speed,
              keyframes[keyframe].delta_time);
        }

        if (keyframes.empty()) camera = config;  // back to start

        if (hasCtrl) {
          // Start playing: spline and start at keyframe.
          if (!keyframes.empty()) {
            CatmullRom(keyframes, &splines, config.loop);
            int nkeys = keyframe;
            for (splines_index = 0; splines_index < splines.size();
                 ++splines_index) {
              if (splines[splines_index].isKey())
              if (nkeys-- == 0) break;
            }
            render_time = splines[splines_index].time;
            render_start = now() - render_time;
          }
        } else {
          splines.clear();
        }
      } break;

      case SDLK_BACKSPACE: {
        --keyframe;
        if (keyframe >= keyframes.size()) keyframe = keyframes.size() - 1;
        if (keyframe < keyframes.size()) {
          camera = keyframes[keyframe];
          printf("at keyframe %lu, speed %.8e, delta_time %f\n",
              (unsigned long)keyframe, keyframes[keyframe].speed,
              keyframes[keyframe].delta_time);
        } else camera = config;
        splines.clear();
      } break;

      // Resolve controller value changes that happened during rendering.
      case SDLK_LEFT:  ctlXChanged = 1; updateControllerX(ctl, -(consecutiveChanges=1), hasAlt); break;
      case SDLK_RIGHT: ctlXChanged = 1; updateControllerX(ctl,  (consecutiveChanges=1), hasAlt); break;
      case SDLK_DOWN:  ctlYChanged = 1; updateControllerY(ctl, -(consecutiveChanges=1), hasAlt); break;
      case SDLK_UP:    ctlYChanged = 1; updateControllerY(ctl,  (consecutiveChanges=1), hasAlt); break;

        // Current keyframe manouvering in screenspace.
        case SDLK_KP4: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].moveAbsolute(camera.right(),
                                             -.5*keyframes[keyframe].speed);
          }
        } break;
        case SDLK_KP6: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].moveAbsolute(camera.right(),
                                             .5*keyframes[keyframe].speed);
          }
        } break;
        case SDLK_KP8: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].moveAbsolute(camera.up(),
                                             .5*keyframes[keyframe].speed);
          }
        } break;
        case SDLK_KP2: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].moveAbsolute(camera.up(),
                                             -.5*keyframes[keyframe].speed);
          }
        } break;

        // Adjust keyframe speed / eye separation. Like mousewheel.
        case SDLK_KP9: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].speed *= 1.1;
          }
        } break;
        case SDLK_KP7: {
          if (keyframe < keyframes.size()) {
            keyframes[keyframe].speed *= .9;
          }
        } break;

        // See whether the active controller has changed.
        default: {
          Controller oldCtl = ctl;
          changeController(event.key.keysym.sym, &ctl);
          if (ctl != oldCtl) { consecutiveChanges = 0; }
        } break;
      }}
      break;
     }
    }

    if (done) break;

    // Get keyboard and mouse state.
    Uint8* keystate = SDL_GetKeyState(0);
    int mouse_dx, mouse_dy;
    Uint8 mouse_buttons = SDL_GetRelativeMouseState(&mouse_dx, &mouse_dy);
    int mouse_button_left = mouse_buttons & SDL_BUTTON(SDL_BUTTON_LEFT);
    int mouse_button_right = mouse_buttons & SDL_BUTTON(SDL_BUTTON_RIGHT);

    bool hasAlt = keystate[SDLK_RALT] || keystate[SDLK_LALT];
    bool hasCtrl = keystate[SDLK_RCTRL] || keystate[SDLK_LCTRL];
    (void)hasCtrl;

    // Continue after calling SDL_GetRelativeMouseState() so view direction
    // does not jump after closing AntTweakBar.
    if (!grabbedInput) continue;

    (void)mouse_buttons;
    (void)mouse_button_left;
    (void)mouse_button_right;

    if (keystate[SDLK_w]) camera.move(0, 0,  camera.speed);  //forward
    if (keystate[SDLK_s]) camera.move(0, 0, -camera.speed);  //back

    if (keystate[SDLK_a]) camera.move(-camera.speed, 0, 0);  //left
    if (keystate[SDLK_d]) camera.move( camera.speed, 0, 0);  //right

    // Mouse look.
    if (grabbedInput && (mouse_dx != 0 || mouse_dy != 0)) {
      m_rotateX2(camera.mouse_rot_speed * mouse_dx * camera.fov_x / 90.0);
      m_rotateY2(camera.mouse_rot_speed * mouse_dy * camera.fov_y / 75.0);
    }
    if (keystate[SDLK_q]) m_rotateZ2(camera.keyb_rot_speed);
    if (keystate[SDLK_e]) m_rotateZ2(-camera.keyb_rot_speed);

   if (keystate[SDLK_z]){ if (camera.speed > 0.000001) camera.speed -= camera.speed/10; printf("speed %.8e\n", camera.speed);}
   if (keystate[SDLK_c]){ if (camera.speed < 1.0) camera.speed += camera.speed/10; printf("speed %.8e\n", camera.speed);}

    // Change the value of the active controller.
    if (!ctlXChanged) {
      if (keystate[SDLK_LEFT])  { ctlXChanged = 1; updateControllerX(ctl, -++consecutiveChanges, hasAlt); }
      if (keystate[SDLK_RIGHT]) { ctlXChanged = 1; updateControllerX(ctl,  ++consecutiveChanges, hasAlt); }
    }
    if (!ctlYChanged) {
      if (keystate[SDLK_DOWN])  { ctlYChanged = 1; updateControllerY(ctl, -++consecutiveChanges, hasAlt); }
      if (keystate[SDLK_UP])    { ctlYChanged = 1; updateControllerY(ctl,  ++consecutiveChanges, hasAlt); }
    }

    if (!(ctlXChanged || ctlYChanged)) consecutiveChanges = 0;
  }

  TwTerminate();

  if (!rendering && keyframesChanged) {
    // TODO: ask whether to save keyframes
  }
  camera.saveConfig("last.cfg", &defines);  // Save a config file on exit, just in case.
  return 0;
}
