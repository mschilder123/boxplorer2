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

#define _strdup strdup
#define __FUNCTION__ "boxplorer"
#define MAX_PATH 256

#else  // _WIN32

#include <winsock2.h>
#include <Ws2tcpip.h>
#include <CommDlg.h>
#include <direct.h>

#pragma warning(disable: 4996) // unsafe function
#pragma warning(disable: 4244) // double / float conversion
#pragma warning(disable: 4305) // double / float truncation
#pragma warning(disable: 4800) // forcing value to bool

#pragma comment(lib, "Ws2_32.lib")
#pragma comment(lib, "SDL.lib")
#pragma comment(lib, "SDLmain.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "shell32.lib")
#pragma comment(lib, "comdlg32.lib")

#if defined(HYDRA)
#include <sixense.h>
#pragma comment(lib, "sixense.lib")
#pragma comment(lib, "sixense_utils.lib")
#endif

typedef SOCKET socket_t;

#endif  // _WIN32

#include <vector>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

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

#include "interpolate.h"
#include "uniforms.h"
#include "camera.h"
#include "shader.h"

#define DEFAULT_CONFIG_FILE "default.cfg"
#define DEFAULT_CONFIG  "cfgs/rrrola/" DEFAULT_CONFIG_FILE
#define VERTEX_SHADER_FILE   "vertex.glsl"
#define FRAGMENT_SHADER_FILE "fragment.glsl"
#define EFFECTS_VERTEX_SHADER_FILE   "effects_vertex.glsl"
#define EFFECTS_FRAGMENT_SHADER_FILE "effects_fragment.glsl"

#define die(...)    ( fprintf(stderr, __VA_ARGS__), exit(-1), 1 )
#ifndef ARRAYSIZE
#define ARRAYSIZE(x) ( sizeof(x)/sizeof((x)[0]) )
#endif

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

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

}  // namespace GLSL


#pragma warning(disable: 4244) // conversion loss
#pragma warning(disable: 4305) // truncation

#define sign(a) GLSL::sign(a)

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

// Our SDL handle.
SDL_Surface* screen;

// Optional #defines for glsl compilation from .cfg file.
string defines;

// Pinhole camera modes.
enum StereoMode { ST_NONE=0,
                  ST_OVERUNDER,
                  ST_XEYED,
                  ST_INTERLACED,
                  ST_SIDEBYSIDE,
                  ST_QUADBUFFER,
                  ST_OCULUS
} stereoMode = ST_NONE;

// ogl framebuffer object, one for each eye.
GLuint fbo[1];
// texture that frame got rendered to
GLuint texture[1];
// depth buffer attached to fbo
GLuint depthBuffer[1];

string BaseDir;     // Where our executable and include dirs live.
string WorkingDir;  // Where current fractal code & data lives.
string BaseFile;     // Initial argument filename.

Shader fractal;
Shader effects;

Uniforms uniforms;

// texture holding background image
GLuint background_texture;

// DLP-Link L/R or interlaced polarity
int polarity=1;

int kJOYSTICK;
SDL_Joystick *joystick;

////////////////////////////////////////////////////////////////
// Helper functions


// Allocate a char[] and read a text file into it. Return 0 on error.
char* _readFile(char const* name) {
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

bool readFile(const string& name, string* content) {
  string filename(WorkingDir + name);
  char* s = _readFile(filename.c_str());
  if (!s) return false;
  content->assign(s);
  free(s);
  printf(__FUNCTION__ " : read '%s'\n", filename.c_str());
  return true;
}


////////////////////////////////////////////////////////////////
// FPS tracking.

int framesToAverage;
Uint32* frameDurations;
int frameDurationsIndex = 0;
Uint32 lastFrameTime;

double now() {
  return (double)SDL_GetTicks() / 1000.0;
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
  float fps = framesToAverage * 1000.f / sum;
#if 1
  static Uint32 lastfps = 0;
  if (lastFrameTime - lastfps > 1000) {
    printf("fps %f\n", fps);
    lastfps = lastFrameTime;
  }
#endif
  return fps;
}

////////////////////////////////////////////////////////////////
// Current logical state of the program.

#include "params.h"

int config_width;
int config_height;

class Camera : public KeyFrame {
  private:


  public:
   Camera& operator=(const KeyFrame& other) {
    *((KeyFrame*)this) = other;
    return *this;
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
      double z_near = fabs(speed);
      double z_far = speed * 65535.0;
      double fH = tan( fov_y * PI / 360.0f ) * z_near;
      double fW = tan( fov_x * PI / 360.0f ) * z_near;
      glFrustum(-fW, fW, -fH, fH, z_near, z_far);

      orthogonalize();
      double matrix[16] = {
         right()[0], up()[0], -ahead()[0], 0,
         right()[1], up()[1], -ahead()[1], 0,
         right()[2], up()[2], -ahead()[2], 0,
                  0,       0,           0, 1
      };
      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixd(matrix);
      // Do not translate, keep eye at 0 to retain drawing precision.
      //glTranslated(-pos()[0], -pos()[1], -pos()[2]);
   }

   // Load configuration.
   bool loadConfig(const string& configFile, string* defines = NULL) {
     bool result = false;
   string filename(WorkingDir + configFile);
     FILE* f;
     if ((f = fopen(filename.c_str(), "r")) != 0) {
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

     // Parse common parameters.
#define PROCESS(type, name, nameString, doSpline) \
           if (!strcmp(s, nameString)) { v=fscanf(f, " %lf", &val); name = val; continue; }
         PROCESS_COMMON_PARAMS
#undef PROCESS

         for (i=0; i<ARRAYSIZE(par); i++) {
           char p[256];
           sprintf(p, "par%lu", (unsigned long)i);
           if (!strcmp(s, p)) {
             v=fscanf(f, " %f %f %f", &par[i][0], &par[i][1], &par[i][2]);
             break;
           }
         }
       }
       fclose(f);
       printf(__FUNCTION__ " : read '%s'\n", configFile.c_str());
       result = true;
     } else {
       printf(__FUNCTION__ " : failed to open '%s'\n", configFile.c_str());
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

     orthogonalize();
   mat2quat(this->v, this->q);

     // Don't do anything with user parameters - they must be
     // sanitized (clamped, ...) in the shader.
   }

   // Save configuration.
   void saveConfig(const string& configFile, string* defines = NULL) {
     FILE* f;
   string filename(WorkingDir + configFile);
     if ((f = fopen(filename.c_str(), "w")) != 0) {
       if (defines != NULL)
         fprintf(f, "%s", defines->c_str());

     // Write common parameters.
#define PROCESS(type, name, nameString, doSpline) \
     fprintf(f, nameString " %g\n", (double)name);
       PROCESS_COMMON_PARAMS
#undef PROCESS

       fprintf(f, "position %12.12e %12.12e %12.12e\n", pos()[0], pos()[1], pos()[2]);
       fprintf(f, "direction %g %g %g\n", ahead()[0], ahead()[1], ahead()[2]);
       fprintf(f, "upDirection %g %g %g\n", up()[0], up()[1], up()[2]);
       for (size_t i=0; i<ARRAYSIZE(par); i++) {
         fprintf(f, "par%lu %g %g %g\n", (unsigned long)i, par[i][0], par[i][1], par[i][2]);
       }
       fclose(f);
       printf(__FUNCTION__ " : wrote '%s'\n", filename.c_str());
     }
   }

   // Send parameters to gpu.
   void setUniforms(float x_scale, float x_offset,
                    float y_scale, float y_offset,
                    double spd) {
     #define glSetUniformf(name) \
       glUniform1f(glGetUniformLocation(program, #name), name);
     #define glSetUniformfv(name) \
       glUniform3fv(glGetUniformLocation(program, #name), ARRAYSIZE(name), (float*)name);
     #define glSetUniformi(name) \
       glUniform1i(glGetUniformLocation(program, #name), name);

   GLuint program = fractal.program();

   // These might be dupes w/ uniforms.send() below.
   // Leave for now until all .cfg got updated.
     glSetUniformi(max_steps); glSetUniformf(min_dist);
     glSetUniformi(iters); glSetUniformi(color_iters);
     glSetUniformf(ao_eps); glSetUniformf(ao_strength);
     glSetUniformf(glow_strength); glSetUniformf(dist_to_color);
   glSetUniformi(nrays); glSetUniformf(focus);

   // Non-user uniforms.
     glSetUniformf(fov_x); glSetUniformf(fov_y);

     glSetUniformf(x_scale); glSetUniformf(x_offset);
     glSetUniformf(y_scale); glSetUniformf(y_offset);

     glSetUniformf(time);

     glUniform1f(glGetUniformLocation(program, "speed"), spd);
     glUniform1f(glGetUniformLocation(program, "xres"), config_width);
     glUniform1f(glGetUniformLocation(program, "yres"), config_height);

     #if defined(GL_ARB_gpu_shader_fp64)
     // Also pass in some double precision values, if supported.
     glUniform1d(glGetUniformLocation(program, "dspeed"), spd);
     glUniform3dv(glGetUniformLocation(program, "deye"), 3, pos());
     #endif

   // Old-style par[] list.
     glSetUniformfv(par);

     #undef glSetUniformf
     #undef glSetUniformfv
     #undef glUniform1i

   // New-style discovered && active uniforms only.
   uniforms.send(program);
   }

   void render(enum StereoMode stereo) {
     activate();  // Load view matrix for shader.
     switch(stereo) {
       case ST_OVERUNDER: {  // left / right
         setUniforms(1.0, 0.0, 2.0, 1.0, +speed);
         glRects(-1,-1,1,0);  // draw bottom half of screen
         setUniforms(1.0, 0.0, 2.0, -1.0, -speed);
         glRects(-1,0,1,1);  // draw top half of screen
         } break;
       case ST_QUADBUFFER: {  // left - right
         glDrawBuffer(GL_BACK_LEFT);
         setUniforms(1.0, 0.0, 1.0, 0.0, -speed*polarity);
         glRects(-1,-1,1,1);
         glDrawBuffer(GL_BACK_RIGHT);
         setUniforms(1.0, 0.0, 1.0, 0.0, +speed*polarity);
         glRects(-1,-1,1,1);
         } break;
       case ST_XEYED: {  // right | left
         setUniforms(2.0, +1.0, 1.0, 0.0, +speed);
         glRectf(-1,-1,0,1);  // draw left half of screen
         setUniforms(2.0, -1.0, 1.0, 0.0, -speed);
         glRectf(0,-1,1,1);  // draw right half of screen
         } break;
       case ST_SIDEBYSIDE: {  // left | right
         setUniforms(2.0, +1.0, 1.0, 0.0, -speed);
         glRectf(-1,-1,0,1);  // draw left half of screen
         setUniforms(2.0, -1.0, 1.0, 0.0, +speed);
         glRectf(0,-1,1,1);  // draw right half of screen
         } break;
       case ST_NONE:
         setUniforms(1.0, 0.0, 1.0, 0.0, speed);
         // Draw screen in N (==2 for now) steps for x-fire / sli
         glRects(-1,-1,0,1);  // draw left half
         glRects(0,-1,1,1);  // draw right half
         break;
       case ST_INTERLACED:
         setUniforms(1.0, 0.0, 1.0, 0.0, speed*polarity);
         glRects(-1,-1,0,1);  // draw left half
         glRects(0,-1,1,1);  // draw right half
         break;
       case ST_OCULUS:
         setUniforms(2.0, +1.0, 1.0, 0.0, -speed);
         glRectf(-1,-1,0,1);  // draw left half of screen
         setUniforms(2.0, -1.0, 1.0, 0.0, +speed);
         glRectf(0,-1,1,1);  // draw right half of screen
         break;
      }
   }

#if defined(HYDRA)
  void mixHydraOrientation(float* quat) {
    double q[4];
    q[0] = quat[0];
    q[1] = quat[1];
    q[2] = quat[2];
    q[3] = quat[3];
    qnormalize(q);
    qmul(q, this->q);
    quat2mat(q, this->v);
  }
  void mixSensorOrientation(socket_t sock) {
    sockaddr_in SenderAddr;
    int SenderAddrSize = sizeof (SenderAddr);

    double q1[4];
  bool gotData = false;

  for(;;) {
      float qf[4];  // receive Quatf
      int r = recvfrom(sock, (char*)&qf, sizeof(qf), 0, (SOCKADDR*)&SenderAddr, &SenderAddrSize);
      if (r == sizeof(qf)) {
        for (int i = 0; i < 4; ++i) q1[i] = qf[i];
    gotData = true;
    } else break;
  }

  if (gotData) {
        q1[2] = -q1[2];  // We roll other way
    qnormalize(q1);

    // combine current view quat with sensor quat.
    qmul(q1, this->q);
        quat2mat(q1, this->v);
  }
  }
#endif

} camera,  // Active camera view.
  config;  // Global configuration set.

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

#define NSUBFRAMES 100  // # splined subframes between keyframes.

void CatmullRom(const vector<KeyFrame>& keyframes,
                vector<KeyFrame>* output,
                bool loop = false,
                int nsubframes = NSUBFRAMES) {
  output->clear();
  if (keyframes.size() < 2) return;  // Need at least two points.

  vector<KeyFrame> controlpoints(keyframes);

  size_t n = controlpoints.size();

  // Compute / check quats.
  // Pick smaller angle between two quats.
  mat2quat(controlpoints[0].v, controlpoints[0].q);
  quat2x(controlpoints[0].q, controlpoints[0].x);
  for (size_t i = 1; i < n; ++i) {
    mat2quat(controlpoints[i].v, controlpoints[i].q);
    double dot =
      controlpoints[i - 1].q[0] * controlpoints[i].q[0] +
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
    if (controlpoints[controlpoints.size()-1].delta_time == 0) {
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
  controlpoints[0].time = 0;   // time starts at 0.
  for (size_t i = 1; i < n - 1; ++i) {
    time += controlpoints[i].delta_time;
    controlpoints[i].time = time;
  }
  // Last one's time is p0 for spline of first one; set to 0 as well.
  controlpoints[n - 1].time = 0;

  // Now spline all into intermediate frames.
  for (size_t i = 0; i < n - 3; ++i) {
    const KeyFrame *p0 = &controlpoints[i>0?i-1:n-1];
    const KeyFrame *p1 = &controlpoints[i];
    const KeyFrame *p2 = &controlpoints[i+1];
    const KeyFrame *p3 = &controlpoints[i+2];
    for (int f = 0; f < nsubframes; ++f) {
      KeyFrame tmp = config;  // Start with default values.
      tmp.setKey(f == 0);
      const double t = ((double)f) / nsubframes;

      // The CatmullRom spline function; 0 <= t <= 1
      // Suffers from overshoot for non-evenly spaced control points.
      // Need to look into Bessel-Overhauser mitigation.
      #define SPLINE(X,p0,p1,p2,p3) \
        ((X) = (double)(.5 * (2 * (p1) + \
                            t*( (-(p0) + (p2)) + \
                                t*( (2*(p0) - 5*(p1) + 4*(p2) - (p3)) + \
                                    t*(-(p0) + 3*(p1) - 3*(p2) + (p3)) ) ) ) ) )

      // Spline over splinable representation of quat.
      for (size_t j = 0; j < 4; ++j) {
        SPLINE(tmp.x[j], p0->x[j], p1->x[j], p2->x[j], p3->x[j]);
      }
      x2quat(tmp.x, tmp.q);  // convert back to quat
    qnormalize(tmp.q);
      quat2mat(tmp.q, tmp.v);  // convert quat to the splined rotation matrix

      // Spline position into tmp.v[12..15]
      for (size_t j = 12; j < 15; ++j) {
        // To control numerical precision, re-base to (p2-p1)/2.
        double a = p0->v[j], b = p1->v[j], c = p2->v[j], d = p3->v[j];
        double base = .5 * (c - b);
        a -= base; b -= base; c -= base; d -= base;
        SPLINE(tmp.v[j], a, b, c, d);
        tmp.v[j] += base;
      }

      // Spline par[] array. Some of those could also be rotations, which will not
      // spline nicely at all times..
      // TODO: have couple of uniform quats for shader use and spline those nicely.
      for (size_t j = 0; j < ARRAYSIZE(tmp.par); ++j) {
        SPLINE(tmp.par[j][0],
               p0->par[j][0], p1->par[j][0], p2->par[j][0], p3->par[j][0]);
        SPLINE(tmp.par[j][1],
               p0->par[j][1], p1->par[j][1], p2->par[j][1], p3->par[j][1]);
        SPLINE(tmp.par[j][2],
               p0->par[j][2], p1->par[j][2], p2->par[j][2], p3->par[j][2]);
      }

    // Spline generic float uniform array.
    for (int j = 0; j < tmp.n_funis; ++j) {
      SPLINE(tmp.funis[j], p0->funis[j], p1->funis[j], p2->funis[j], p3->funis[j]);
    }

      // Spline common params, if marked as such.
#define PROCESS(a,b,c,doSpline) \
      if (doSpline) { SPLINE(tmp.b, p0->b, p1->b, p2->b, p3->b); }
      PROCESS_COMMON_PARAMS
#undef PROCESS

      #undef SPLINE

      //tmp.orthogonalize();  // this should be no-op given matrix came from quat?
      output->push_back(tmp);
    }
  }
}

////////////////////////////////////////////////////////////////
// Controllers.

typedef enum Controller {
  // user parameters: 0..9
  // other parameters:
  CTL_FOV = ARRAYSIZE(camera.par), CTL_RAY, CTL_ITER, CTL_AO, CTL_GLOW,
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
void m_progressiveAdd(double* x, int d) { *x += 0.001 * (sign(d) * ((abs(d)+4) / 4)); }
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
        y, camera.par[c][1],
        x, camera.par[c][0]
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
  string vertex(default_vs);
  string fragment(default_fs);

  readFile(VERTEX_SHADER_FILE, &vertex);
  readFile(FRAGMENT_SHADER_FILE, &fragment);

  glsl_source.assign(defines + fragment);

  return fractal.compile(defines, vertex, fragment);
}

// Compile and activate shader programs for frame buffer manipulation.
// Return the program handle.
int setupShaders2(void) {
  string vertex(effects_default_vs);
  string fragment(effects_default_fs);

  readFile(EFFECTS_VERTEX_SHADER_FILE, &vertex);
  readFile(EFFECTS_FRAGMENT_SHADER_FILE, &fragment);

  glsl_source.append(fragment);

  return effects.compile(defines, vertex, fragment);
}

bool setupDirectories(const char* configFile) {
  BaseFile.clear();
  BaseDir.clear();
  WorkingDir.clear();

  char dirName[MAX_PATH];
#if defined(_WIN32)
  GetModuleFileName(NULL, dirName, MAX_PATH);
  BaseDir.assign(dirName);
  size_t slash = BaseDir.rfind('\\');
  if (slash != string::npos) BaseDir.erase(slash + 1);

  char* fileName = NULL;
  DWORD result = GetFullPathName(configFile, MAX_PATH, dirName, &fileName);
  if (result) {
    if (fileName) {
      BaseFile.assign(fileName);
      *fileName = '\0';
    }
    WorkingDir.assign(dirName);
  }
#else
  strncpy(dirName, configFile, sizeof dirName);
  dirName[sizeof dirName - 1] = 0;
  int i = strlen(dirName);
  while (i > 0 && strchr("/", dirName[i - 1]) == NULL) --i;
  BaseFile.assign(&dirName[i]);
  dirName[i] = 0;
  BaseDir.assign("./");  // Assume relative to cwd.
  WorkingDir.assign(dirName);
#endif

  if (BaseFile.empty()) BaseFile.assign(DEFAULT_CONFIG_FILE);

  cout << __FUNCTION__ << " : " << BaseDir
                       << ", " << WorkingDir
                       << ", " << BaseFile << endl;

  return true;
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

  if(stereoMode==ST_QUADBUFFER) {
    SDL_GL_SetAttribute(SDL_GL_STEREO, 1);
  }

  (screen = SDL_SetVideoMode(config.width, config.height, bpp,
      SDL_OPENGL | (config.fullscreen ? SDL_FULLSCREEN : SDL_RESIZABLE)))
    || die("Video mode initialization failed: %s\n", SDL_GetError());

  if(stereoMode==ST_QUADBUFFER) {
    int ga = 0;
    SDL_GL_GetAttribute(SDL_GL_STEREO, &ga);
    if (ga == 0) die("No stereo rendering available: %s\n", SDL_GetError());
  }

  config_width = config.width;
  config_height = config.height;

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

  (setupShaders()) ||
      die("Error in GLSL fractal shader compilation:\n%s\n", fractal.log().c_str());

  if (background.data() != NULL) {
    // Load background image into texture
    glDeleteTextures(1, &background_texture);  // free existing
    glGenTextures(1, &background_texture);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, background_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8,
                 background.width(), background.height(),
                 0, GL_BGR, GL_UNSIGNED_BYTE, background.data());
    glGenerateMipmap(GL_TEXTURE_2D);
    printf(__FUNCTION__ " : background texture at %d\n", background_texture);
    glBindTexture(GL_TEXTURE_2D, 0);
  }

  if ((status = glGetError()) != GL_NO_ERROR)
    die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

  if (config.enable_dof) {
    // Compile DoF shader, setup FBO as render target.

    (setupShaders2()) ||
        die("Error in GLSL effects shader compilation:\n%s\n", effects.log().c_str());

    // Create depth buffer(s)
    glDeleteRenderbuffers(ARRAYSIZE(depthBuffer), depthBuffer);  // free existing
    glGenRenderbuffers(ARRAYSIZE(depthBuffer), depthBuffer);

    // Create textures to render to
    glDeleteTextures(ARRAYSIZE(texture), texture);  // free existing
    glGenTextures(ARRAYSIZE(texture), texture);

    // Create framebuffers
    glDeleteFramebuffers(ARRAYSIZE(fbo), fbo);  // free existing
    glGenFramebuffers(ARRAYSIZE(fbo), fbo);

  for (int i = 0; i < ARRAYSIZE(fbo); ++i) {
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer[i]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                          config.width, config.height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    if ((status = glGetError()) != GL_NO_ERROR)
      die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

    glBindTexture(GL_TEXTURE_2D, texture[i]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);

    // Allocate storage, float rgba if available. Pick max alpha width.
#ifdef GL_RGBA32F
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, config.width, config.height,
                 0, GL_BGRA, GL_FLOAT, NULL);
#else
#ifdef GL_RGBA16
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, config.width, config.height,
                 0, GL_BGRA, GL_UNSIGNED_SHORT, NULL);
#else
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, config.width, config.height,
                 0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);
    fprintf(stderr, __FUNCTION__ " : 8 bit alpha, very granular DoF experience ahead..\n");
#endif
#endif

    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo[i]);

    // Attach texture to framebuffer as colorbuffer
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           texture[i], 0);

    // Attach depthbuffer to framebuffer
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, depthBuffer[i]);

    if ((status = glGetError()) != GL_NO_ERROR)
      die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);

    status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE)
      die(__FUNCTION__ " : glCheckFramebufferStatus() : %04x\n", status);

    // Back to normal framebuffer.
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

    if ((status = glGetError()) != GL_NO_ERROR)
      die(__FUNCTION__ "[%d] : glGetError() : %04x\n", __LINE__, status);
  }
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
    if (index < 0 || index > (int)ARRAYSIZE(camera.par)) continue;
    if (_xyz < 'x' || _xyz > 'z') continue;

    printf("parameter %s par[%d].%c {%s}\n",
           varName.c_str(), index, _xyz, attr.c_str());

    float* address = &camera.par[index][_xyz-'x'];

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

  if (stereoMode == ST_OCULUS) {
  // Position HUD center for left eye.
    char pos[100];
  int x = config.width;
  int y = config.height;
  sprintf(pos, "boxplorer position='%d %d'", x/6, y/4);
    TwDefine(pos);
  TwDefine("GLOBAL fontsize=3");
  }

#if 0
  // Add TW UI for common parameters.
#define PROCESS(a,b,c,d) initTwUniform(c, &camera.b);
  PROCESS_COMMON_PARAMS
#undef PROCESS
#else
  uniforms.bindToUI(bar);
#endif

  initTwParDefines();
}

void LoadKeyFrames(bool fixedFov) {
  char filename[256];
  for (int i = 0; ; ++i) {
    sprintf(filename, "%s-%u.cfg", kKEYFRAME, i);
  // We load into global camera since that's where
  // the uniforms are bound to.
    if (!camera.loadConfig(filename)) break;
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

////////////////////////////////////////////////////////////////
// Setup, input handling and drawing.

int main(int argc, char **argv) {
  bool rendering = false;
  bool loop = false;
  bool useTime = false;
  bool configSpeed = false;
  bool fixedFov = false;
  int enableDof = 0;

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
    } else if (!strcmp(argv[argc-1], "--quadbuffer")) {
      stereoMode = ST_QUADBUFFER;
    } else if (!strcmp(argv[argc-1], "--oculus")) {
      stereoMode = ST_OCULUS;
      defines.append("#define ST_OCULUS\n");
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
    } else if (!strncmp(argv[argc-1], "--joystick=", 11)) {
      kJOYSTICK = atoi(argv[argc-1] + 11);
    } else break;
    --argc;
  }

  if (stereoMode == ST_NONE) defines.append("#define ST_NONE\n");

  const char* configFile = (argc>=2 ? argv[1] : DEFAULT_CONFIG);

  // Load configuration.
  if (setupDirectories(configFile) &&
    config.loadConfig(BaseFile, &defines)) {
  } else {
#if defined(_WIN32)
    // For windows users that don't specify an argument, offer a dialog.
    char filename[256] = {0};
    OPENFILENAME opf = {0};
    opf.lpstrFilter = "Fractal configuration files\0*.cfg\0\0";
    opf.nFilterIndex = 1L;
    opf.lpstrFile = filename;
    opf.nMaxFile = 255;
    opf.nMaxFileTitle=50;
    opf.lpstrTitle = "Open fractal configuration file";
    opf.nFileOffset = 0;
    opf.lpstrDefExt = "cfg";
    opf.Flags = OFN_PATHMUSTEXIST | OFN_READONLY;
    opf.lStructSize = sizeof(OPENFILENAME);
    if (GetOpenFileName(&opf) &&
    setupDirectories(opf.lpstrFile) &&
    config.loadConfig(BaseFile, &defines)) {
    } else
#endif
    die("Usage: boxplorer <configuration-file.cfg>\n");
  }

#if defined(HYDRA)
  // Listen on UDP:1337 for quat sent from Oculus SensorBox
  WSADATA wsaData;
  WSAStartup(MAKEWORD(2, 2), &wsaData);

  socket_t RecvSocket = INVALID_SOCKET;

  if (stereoMode == ST_OCULUS) {
    RecvSocket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    u_long nonblock = 1;
    ioctlsocket(RecvSocket, FIONBIO, &nonblock);
    sockaddr_in RecvAddr;
    RecvAddr.sin_family = AF_INET;
    RecvAddr.sin_port = htons(1337);
    RecvAddr.sin_addr.s_addr = htonl(INADDR_ANY);
    bind(RecvSocket, (SOCKADDR *) & RecvAddr, sizeof (RecvAddr));
  }

  if (sixenseInit() != SIXENSE_SUCCESS) {
    die("sixenseInit() fail!");
  }

  sixenseSetFilterEnabled(1);
  if (sixenseSetFilterParams(0.0, 0.0, 1000.0, 1.0) != SIXENSE_SUCCESS) {
    die("SetFilterParams() fail!");
  }

  if (sixenseSetActiveBase(0) != SIXENSE_SUCCESS) {
    die("sixenseSetActiveBase() fail!");
  }

  sixenseAllControllerData ssdata;

  if (sixenseIsControllerEnabled(0) != SIXENSE_SUCCESS) {
    die("controller(0) not enabled!");
  }
#endif

  // Sanitize / override config parameters.
  if (loop) config.loop = true;
  if (enableDof) config.enable_dof = (enableDof == 1);  // override
  if (stereoMode == ST_INTERLACED || stereoMode == ST_QUADBUFFER) config.enable_dof = 0;  // mipmapping does not work for interlaced.
  if (stereoMode == ST_OCULUS) {
    config.width = 1280; config.height = 800;  // Fix rez. Otherwise mirrored screen drops Rift?
  config.fov_y = 110; config.fov_x = 90.0;
  fixedFov = true;
  }
  if (config.fps < 5) config.fps = 30;
  if (config.depth_size < 16) config.depth_size = 16;
  if (stereoMode == ST_XEYED) config.width *= 2;

  int savedWidth = config.width;
  int savedHeight = config.height;

  LoadBackground();

  // Initialize SDL and OpenGL graphics.
  SDL_Init(SDL_INIT_VIDEO) == 0 ||
      die("SDL initialization failed: %s\n", SDL_GetError());
  atexit(SDL_Quit);

  if(kJOYSTICK) {
    SDL_InitSubSystem(SDL_INIT_JOYSTICK);
    for(int i=0; i < SDL_NumJoysticks(); i++) printf(__FUNCTION__ " : JoystickName %i: '%s'\n", i+1, SDL_JoystickName(i));
    joystick = SDL_JoystickOpen(kJOYSTICK-1);
    printf(__FUNCTION__ " : JoystickNumAxes   : %i\n", SDL_JoystickNumAxes(joystick));
    printf(__FUNCTION__ " : JoystickNumButtons: %i\n", SDL_JoystickNumButtons(joystick));
    printf(__FUNCTION__ " : JoystickNumHats   : %i\n", SDL_JoystickNumHats(joystick));
  }

   // Set up the video mode, OpenGL state, shaders and shader parameters.
  initGraphics();

  // Parse as many uniforms from glsl source as we can find.
  uniforms.parseFromGlsl(glsl_source);

  // TODO: prune uniforms to just those reported active by shader compiler.
  cout << fractal.uniforms();

  // Bind as many uniforms as we can find a match for to camera.
  uniforms.link(&camera);

  // Load key frames, if any. Clobbers camera.
  LoadKeyFrames(fixedFov);

  // Load initial camera; sets all known, linked uniforms.
  camera.loadConfig(BaseFile);

  initTwBar();
  initFPS(FPS_FRAMES_TO_AVERAGE);

  //printf(__FUNCTION__ " : GL_EXTENSIONS: %s\n", glGetString(GL_EXTENSIONS));

  // Main loop.
  Controller ctl = CTL_CAM;  // the default controller is camera rotation
  int consecutiveChanges = 0;

  int done = 0;
  int frameno = 0;

  vector<KeyFrame> splines;
  size_t splines_index = 0;

  double frame_time = 1 / config.fps;
  printf(__FUNCTION__ " : frame time %g\n", frame_time);
  double render_time = 0;
  double render_start = 0;

  double dist_along_spline = 0;
  size_t keyframe = keyframes.size();

  bool ignoreNextMouseUp = false;
  bool dragging = false;

  if (rendering) {
    // Rendering a sequence to disk. Spline the keyframes now.
    CatmullRom(keyframes, &splines, config.loop);
  }

  double last_de = 10.0;

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

  while (!done) {
    int ctlXChanged = 0, ctlYChanged = 0;

    // Splined keyframes playback logic. Messy. Sets up next camera.
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
            double n = now();
            if (n > render_start + camera.time) continue;  // late, skip frame.
            double w = (render_start + camera.time) - n;
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
    } else {
      camera.time = now();
    }

#if defined(HYDRA)
  if (!rendering) {
      // When not rendering a sequence, now mix in orientation (and translation)
      // external sensors might have to add (e.g. Oculus orientation) into the view
      // we are about to render.
      if (stereoMode == ST_OCULUS) {
      camera.mixSensorOrientation(RecvSocket);
    }
    //camera.mixHydraOrientation(ssdata.controllers[0].rot_quat);
  }
#endif

    if (!rendering && (de_func || de_func_64)) {
      // Get a distance estimate, for navigation and eye separation.
      // Setup just the vars needed for DE. For now, iters and par[0..20]
      // TODO: make de() method of camera?
      GLSL::iters = camera.iters;
      for (int i = 0; i < 20; ++i) {
        GLSL::par[i] = GLSL::vec3(camera.par[i]);
      }

      GLSL::dvec3 pos(camera.pos());
    double de = de_func_64?GLSL::abs(de_func_64(pos)):GLSL::abs(de_func(pos));

      if (de != last_de) {
        printf("de=%12.12e\n", de);
        camera.speed = de / 10.0;
        last_de = de;
      }
    }

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_ALWAYS);  // we're writing Z every pixel

    if (background_texture) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, background_texture);
    }

  GLuint program = fractal.program();
    glUseProgram(program);  // the fractal shader

    glUniform1i(glGetUniformLocation(program, "bg_texture"), 0);
    glUniform1i(glGetUniformLocation(program, "use_bg_texture"),
                background_texture);

   if (config.enable_dof && camera.dof_scale != 0) {
      // Render to different back buffer to compute DoF effects later,
      // using alpha as depth channel.
      glBindFramebuffer(GL_FRAMEBUFFER, fbo[frameno&0]);
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    camera.render(stereoMode);  // draw fractal

    glUseProgram(0);

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // If we're rendering some DoF effects, draw from buffer to screen.
    if (config.enable_dof && camera.dof_scale != 0) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glEnable(GL_TEXTURE_2D);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture[frameno&0]);
      glGenerateMipmap(GL_TEXTURE_2D);  // generate mipmaps of our rendered frame.

    GLuint dof_program = effects.program();
      glUseProgram(dof_program);  // Activate our alpha channel DoF shader.

      glUniform1i(glGetUniformLocation(dof_program, "my_texture"), 0);
      glUniform1f(glGetUniformLocation(dof_program, "dof_scale"),
                  camera.dof_scale);
      glUniform1f(glGetUniformLocation(dof_program, "dof_offset"),
                  camera.dof_offset);
      glUniform1f(glGetUniformLocation(dof_program, "speed"),
                  camera.speed);
      #if defined(GL_ARB_gpu_shader_fp64)
      // Also pass in double precision speed, if supported.
      glUniform1d(glGetUniformLocation(dof_program, "dspeed"), camera.speed);
      #endif

      glUniform1f(glGetUniformLocation(dof_program, "xres"), config.width);
      glUniform1f(glGetUniformLocation(dof_program, "yres"), config.height);

      // Ortho projection, entire screen in regular pixel coordinates.
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(0, config.width, config.height, 0, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      // Draw our texture covering entire screen, running the frame shader.
    // Top half
      glBegin(GL_QUADS);
        glTexCoord2f(0,1);
        glVertex2f(0,0);
        glTexCoord2f(0,0.5);
        glVertex2f(0,config.height/2);
        glTexCoord2f(1,0.5);
        glVertex2f(config.width,config.height/2);
        glTexCoord2f(1,1);
        glVertex2f(config.width,0);
      glEnd();

    // Bottom half
      glBegin(GL_QUADS);
        glTexCoord2f(0,0.5);
        glVertex2f(0,config.height/2);
        glTexCoord2f(0,0);
        glVertex2f(0,config.height);
        glTexCoord2f(1,0);
        glVertex2f(config.width,config.height);
        glTexCoord2f(1,0.5);
        glVertex2f(config.width,config.height/2);
      glEnd();

      glUseProgram(0);

      glDisable(GL_TEXTURE_2D);
    }

    // Draw keyframe splined path, if we have 2+ keyframes and not rendering.
    if (!config.no_spline &&
        stereoMode == ST_NONE &&  // TODO: show path in 3d
        keyframes.size() > 1 &&
        splines.empty()) {
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
        // Translate points, eye is at origin.
        glVertex3d(splines[i].pos()[0] - camera.pos()[0],
                   splines[i].pos()[1] - camera.pos()[1],
                   splines[i].pos()[2] - camera.pos()[2]);

        splines[i].move(-2*splines[i].speed, 0, 0);
      }
      glEnd();
      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // left eye
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
        glColor4f(0,0,1 - (i/256.0),1);  // Encode keyframe # in color.
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

    // Process UI events.
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if(event.type==SDL_JOYBUTTONDOWN) {
        // Hack to reach "case SDLK_*" code:
        event.type=SDL_KEYDOWN;
        if(event.jbutton.button) event.key.keysym.sym=SDLK_TAB;
        else event.key.keysym.sym=SDLK_BACKSPACE;
      }
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
              // Peek at framebuffer color for keyframe markers.
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
             SDL_SetCursor(hand_cursor);
             // Drag currently selected keyframe around.
             if (keyframe < keyframes.size()) {
               // Should really be some screenspace conversion..
               // but this works ok for now.
               double fY = 2.0 * tan(camera.fov_y * PI / 360.0f) / config.height;
               double fX = 2.0 * tan(camera.fov_x * PI / 360.0f) / config.width;
               keyframes[keyframe].moveAbsolute(camera.up(),
                   event.motion.yrel*-fY*camera.distanceTo(keyframes[keyframe]));
               keyframes[keyframe].moveAbsolute(camera.right(),
                   event.motion.xrel*fX*camera.distanceTo(keyframes[keyframe]));
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

      // Switch fullscreen mode (drops the whole OpenGL context in Windows).
      case SDLK_RETURN: case SDLK_KP_ENTER: {
        config.fullscreen ^= 1; grabbedInput = 1;
        if (config.fullscreen) {
          savedWidth = config.width; savedHeight = config.height;
          if (stereoMode == ST_INTERLACED || stereoMode == ST_OVERUNDER) {
            // Zalman or acer
            config.height = 1080; config.width = 1920;
          } else if (stereoMode == ST_SIDEBYSIDE || stereoMode == ST_OVERUNDER) {
            // HMZ-T1
            config.height = 720; config.width = 1280;
          } else if (stereoMode == ST_OCULUS) {
            // Oculus Rift
            config.height = 800; config.width = 1280; // WTF? xfire fail on that rez
            //config.height = 720; config.width = 1280;
            config.fov_y = 110.0; config.fov_x = 90.0;
          } else if (stereoMode == ST_NONE) {
            config.height = 1600; config.width = 2560;  // 30"
//            config.height = 1080; config.width = 1920;  // 27"
          }
        } else {
          config.width = savedWidth; config.height = savedHeight;
        }
        initGraphics(); initTwBar();
      } break;

      // Switch DLP-Link L/R polarity
      case SDLK_p: {
        polarity*=-1;
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
            dist_along_spline = 0;
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

      // Currently selected keyframe manouvering in screenspace.
      // 4-6 left-right, 8-2 up-down, 9-1 further-closer, 7-3 speed at frame.
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
      case SDLK_KP9: {
        if (keyframe < keyframes.size()) {
          keyframes[keyframe].moveAbsolute(camera.ahead(),
                                           .5*keyframes[keyframe].speed);
        }
      } break;
      case SDLK_KP1: {
        if (keyframe < keyframes.size()) {
          keyframes[keyframe].moveAbsolute(camera.ahead(),
                                           -.5*keyframes[keyframe].speed);
        }
      } break;
      case SDLK_KP3: {
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
    Sint16 joystick_x=0, joystick_y=0, joystick_z=0, joystick_r=0;
    Uint8 joystick_hat=0;
    Uint8 mouse_buttons = SDL_GetRelativeMouseState(&mouse_dx, &mouse_dy);
    int mouse_button_left = mouse_buttons & SDL_BUTTON(SDL_BUTTON_LEFT);
    int mouse_button_right = mouse_buttons & SDL_BUTTON(SDL_BUTTON_RIGHT);

    bool hasAlt = keystate[SDLK_RALT] || keystate[SDLK_LALT];
    bool hasCtrl = keystate[SDLK_RCTRL] || keystate[SDLK_LCTRL];
    (void)hasCtrl;

    // Continue after calling SDL_GetRelativeMouseState() so view direction
    // does not jump after closing AntTweakBar.
    if (!grabbedInput) continue;

    if(joystick) {
      SDL_JoystickUpdate();
      joystick_x   = SDL_JoystickGetAxis(joystick, 0);
      joystick_y   = SDL_JoystickGetAxis(joystick, 1);
      joystick_z   = SDL_JoystickGetAxis(joystick, 2);
      joystick_r   = SDL_JoystickGetAxis(joystick, 3);
      joystick_hat = SDL_JoystickGetHat (joystick, 0);
    }

    (void)mouse_buttons;
    (void)mouse_button_left;
    (void)mouse_button_right;

    if (keystate[SDLK_w]) camera.move(0, 0,  camera.speed);  //forward
    if (keystate[SDLK_s]) camera.move(0, 0, -camera.speed);  //back
    if (joystick_z != 0)  camera.move(0, 0,  camera.speed * -joystick_z / 10000.0);

    if (keystate[SDLK_a] || (joystick_hat & SDL_HAT_LEFT )) camera.move(-camera.speed, 0, 0);  //left
    if (keystate[SDLK_d] || (joystick_hat & SDL_HAT_RIGHT)) camera.move( camera.speed, 0, 0);  //right
    if (joystick_hat & SDL_HAT_DOWN) camera.move(0, -camera.speed, 0);  //down
    if (joystick_hat & SDL_HAT_UP  ) camera.move(0,  camera.speed, 0);  //up

    // Mouse look.
    if (grabbedInput && (mouse_dx != 0 || mouse_dy != 0)) {
      m_rotateX2(camera.mouse_rot_speed * mouse_dx * camera.fov_x / 90.0);
      m_rotateY2(camera.mouse_rot_speed * mouse_dy * camera.fov_y / 75.0);
    }
    if (keystate[SDLK_q]) m_rotateZ2(camera.keyb_rot_speed);
    if (keystate[SDLK_e]) m_rotateZ2(-camera.keyb_rot_speed);

    // Joystick look.
    if (joystick_x != 0 || joystick_y != 0) {
      m_rotateX2(camera.mouse_rot_speed *  joystick_x * camera.fov_x / 90.0 / 5000.0);
      m_rotateY2(camera.mouse_rot_speed * -joystick_y * camera.fov_y / 75.0 / 5000.0);
    }
    if (joystick_r != 0) m_rotateZ2(camera.keyb_rot_speed * joystick_r / 10000.0);

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

#if defined(HYDRA)
  // Sixense Hydra
    if (sixenseGetAllNewestData(&ssdata) != SIXENSE_SUCCESS) {
    die("sixenseGetAllNewestData() fail!");
    }
//  printf("%f %f %f %f\n", ssdata.controllers[0].rot_quat[0],ssdata.controllers[0].rot_quat[1],ssdata.controllers[0].rot_quat[2],ssdata.controllers[0].rot_quat[3]);

  camera.move(0, 0, camera.speed *   ssdata.controllers[0].joystick_y);
  m_rotateX2(camera.keyb_rot_speed *.1 * ssdata.controllers[1].joystick_x);
  m_rotateY2(camera.keyb_rot_speed *.1 * ssdata.controllers[1].joystick_y);
  m_rotateZ2(camera.keyb_rot_speed *.1 * -ssdata.controllers[0].joystick_x);

//  printf("%08x, %f\n", ssdata.controllers[0].buttons, ssdata.controllers[0].trigger);
//  printf("%f %f %f\n", ssdata.controllers[0].pos[0],ssdata.controllers[0].pos[1],ssdata.controllers[0].pos[2]);
#endif

    if (!(ctlXChanged || ctlYChanged)) consecutiveChanges = 0;
  }

  TwTerminate();

  camera.width = savedWidth; camera.height = savedHeight;
  camera.saveConfig("last.cfg", &defines);  // Save a config file on exit, just in case.

#if defined(HYDRA)
  closesocket(RecvSocket);
  WSACleanup();

  sixenseExit();
#endif

  return 0;
}
