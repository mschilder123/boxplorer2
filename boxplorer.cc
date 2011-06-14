#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>

#if !defined(_WIN32)

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

using namespace std;

#define NO_SDL_GLEXT
#include <SDL/SDL_opengl.h>
#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#include <SDL/SDL_main.h>

#include <AntTweakBar.h>

#include "shader_procs.h"
#include "default_shaders.h"

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
#define sign(x)     ( (x)<0 ? -1 : 1 )

#define zNear 0.0001f
#define zFar  5.0f

#define FPS_FRAMES_TO_AVERAGE 20

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

// Pinhole camera modes.
enum StereoMode { ST_NONE=0, ST_OVERUNDER, ST_XEYED, ST_INTERLACED }
    stereoMode = ST_NONE;

// ogl framebuffer object
GLuint fbo;
// texture that frame got rendered to
GLuint texture;
// depth buffer attached to fbo
GLuint depthBuffer;
// the dof mipmapper program
int dof_program;

////////////////////////////////////////////////////////////////
// Helper functions

// Compute the dot product of two vectors.
float dot(const float x[3], const float y[3]) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// Normalize a vector. If it was zero, return 0.
int normalize(float x[3]) {
  float len = dot(x, x); if (len == 0) return 0;
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
  PROCESS(float, speed, "speed", true) \
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
  PROCESS(float, asymmetry, "asymmetry", true)

char* parName[10][3];

class KeyFrame {
  public:
   // View matrix.
   float v[16];

   // Declare fractal and other parameters.
   #define PROCESS(a,b,c,d) a b;
   PROCESS_CONFIG_PARAMS
   #undef PROCESS

   // Par[] parameter array.
   float par[20][3];  // min(this, glsl) gets sent to shader.

   bool isKey_;  // Whether this frame is actually a KeyFrame.

   KeyFrame() { memset(this, 0, sizeof *this); }

   float distanceTo(const KeyFrame& other) const {
      float delta[3] = { v[12]-other.v[12],
                         v[13]-other.v[13],
                         v[14]-other.v[14] };
      return sqrt(dot(delta, delta));
   }
   float* right() { return &v[0]; }
   float* up() { return &v[4]; }
   float* ahead() { return &v[8]; }
   float* pos() { return &v[12]; }

   void setKey(bool key) { isKey_ = key; }
   bool isKey() const { return isKey_; }

   void orthogonalize() {
      if (!normalize(ahead())) { ahead()[0]=ahead()[1]=0; ahead()[2]=1; }
      float l = dot(ahead(), up());
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
   void move(float x, float y, float z) {
      for (int i=0; i<3; i++) {
         pos()[i] += right()[i]*x + up()[i]*y + ahead()[i]*z;
      }
   }

   // Move camera in the normalized absolute direction `dir` by `len` units.
   void moveAbsolute(float* dir, float len) {
      for (int i=0; i<3; i++) {
         pos()[i] += len * dir[i];
      }
   }

   // Rotate the camera by `deg` degrees around a normalized axis.
   // Behaves like `glRotate` without normalizing the axis.
   void rotate(float deg, float x, float y, float z) {
     float s = sin(deg*PI/180), c = cos(deg*PI/180), t = 1-c;
     float r[3][3] = {
      { x*x*t +   c, x*y*t + z*s, x*z*t - y*s },
      { y*x*t - z*s, y*y*t +   c, y*z*t + x*s },
      { z*x*t + y*s, z*y*t - x*s, z*z*t +   c }
     };
     for (int i=0; i<3; i++) {
      float c[3];
      for (int j=0; j<3; j++) c[j] = v[i+j*4];
      for (int j=0; j<3; j++) v[i+j*4] = dot(c, r[j]);
     }
   }

   // Set the OpenGL modelview matrix to the camera matrix, for shader.
   void activate() const {
      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixf(v);
   }

   // Set the OpenGL modelview and projection for gl*() functions.
   void activateGl() {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      float fH = tan( fov_y * PI / 360.0f ) * zNear;
      float fW = tan( fov_x * PI / 360.0f ) * zNear;
      glFrustum(-fW, fW, -fH, fH, zNear, zFar);

      orthogonalize();
      float matrix[16] = {
         right()[0], up()[0], -ahead()[0], 0,
         right()[1], up()[1], -ahead()[1], 0,
         right()[2], up()[2], -ahead()[2], 0,
                  0,       0,           0, 1
      };
      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixf(matrix);
      glTranslatef(-pos()[0], -pos()[1], -pos()[2]);
   }

   // Load configuration.
   bool loadConfig(char const* configFile) {
     bool result = false;
     FILE* f;
     if ((f = fopen(configFile, "r")) != 0) {
      size_t i;
      char s[32768];  // max line length
      while (fscanf(f, " %s", s) == 1) {  // read word
        if (s[0] == 0 || s[0] == '#') continue;

        double val;
        int v;

        if (!strcmp(s, "position")) { v=fscanf(f, " %f %f %f", &pos()[0], &pos()[1], &pos()[2]); continue; }
        if (!strcmp(s, "direction")) { v=fscanf(f, " %f %f %f", &ahead()[0], &ahead()[1], &ahead()[2]); continue; }
        if (!strcmp(s, "upDirection")) { v=fscanf(f, " %f %f %f", &up()[0], &up()[1], &up()[2]); continue; }

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
   void saveConfig(char const* configFile) {
     FILE* f;
     if ((f = fopen(configFile, "w")) != 0) {
        #define PROCESS(type, name, nameString, doSpline) \
          fprintf(f, nameString " %g\n", (double)name);
        PROCESS_CONFIG_PARAMS
        #undef PROCESS

        fprintf(f, "position %g %g %g\n", pos()[0], pos()[1], pos()[2]);
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
                    float speed = 0.0) {
     #define glSetUniformf(name) \
       glUniform1f(glGetUniformLocation(program, #name), name);
     #define glSetUniformv(name) \
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
     glSetUniformf(speed);

     glSetUniformv(par);

     #undef glSetUniformf
     #undef glSetUniformv
     #undef glUnifor1i
   }

   void render(enum StereoMode stereo) {
     float p[3] = { pos()[0], pos()[1], pos()[2] };
     switch(stereo) {
       case ST_OVERUNDER:
         move(speed, 0, 0);  // step right
         activate();
         setUniforms(1.0-asymmetry/2, -asymmetry/2, 2.0, 1.0);
         glRects(-1,-1,1,0);  // draw bottom half of screen
         pos()[0] = p[0]; pos()[1] = p[1]; pos()[2] = p[2];  // restore pos
         move(-speed, 0, 0);  // step left
         activate();
         setUniforms(1.0-asymmetry/2, asymmetry/2, 2.0, -1.0);
         glRects(-1,0,1,1);  // draw top half of screen
         pos()[0] = p[0]; pos()[1] = p[1]; pos()[2] = p[2];
         break;
       case ST_XEYED:
         move(speed, 0, 0);  // step right
         activate();
         setUniforms(2.0*(1.0-asymmetry/2), 1.0-asymmetry, 1.0, 0.0);
         glRectf(-1,-1,0,1);  // draw left half of screen
         pos()[0] = p[0]; pos()[1] = p[1]; pos()[2] = p[2];
         move(-speed, 0, 0);  // step left
         activate();
         setUniforms(2.0*(1.0-asymmetry/2), -1.0+asymmetry, 1.0, 0.0);
         glRectf(0,-1,1,1);  // draw right half of screen
         pos()[0] = p[0]; pos()[1] = p[1]; pos()[2] = p[2];
         break;
       case ST_NONE:
         activate();
         setUniforms(1.0, 0.0, 1.0, 0.0);
         glRects(-1,-1,1,1);  // draw entire screen
         break;
       case ST_INTERLACED:
         activate();
         //TODO: add asymmetric frustum to interlaced rendering.
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
    float dist = camera.distanceTo(keyframes[keyframes.size() - 1]);
    float steps = dist / camera.speed;
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
      const float t = f * (1. / nsubframes);

      // The CatmullRom spline function; 0 <= t <= 1
      #define SPLINE(X,p0,p1,p2,p3) \
        (X = (float)(.5 * ( (2 * p1 + \
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
    case CTL_AO: sprintf(s, "aO %.2e d %.2e", camera.ao_strength, camera.ao_eps); break;
    case CTL_GLOW: sprintf(s, "Glow %.3f bgd %.2e", camera.glow_strength, camera.dist_to_color); break;
    case CTL_CAM: {
      sprintf(s, "Look [%4d %4d %4d]",
              (int)(camera.ahead()[0]*100),
              (int)(camera.ahead()[1]*100),
              (int)(camera.ahead()[2]*100));
    } break;
    case CTL_TIME: {
      sprintf(s, "Speed %f DeltaT %.3f", camera.speed, camera.delta_time);
    } break;
    case CTL_3D: {
      sprintf(s, "Sep %f Asym %.3f", camera.speed, camera.asymmetry);
    } break;
  }
  return s;
}

// Update controller.y by the signed count of consecutive changes.
void updateControllerY(Controller c, int d) {
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
    case CTL_3D: m_progressiveAdd(&camera.asymmetry, d); break;
  }
  // Enforce sane bounds.
  if (camera.delta_time < 0) camera.delta_time = 0;
  if (camera.asymmetry < 0) camera.asymmetry = 0;
}

// Update controller.x by the signed count of consecutive changes.
void updateControllerX(Controller c, int d) {
  assert(c <= CTL_LAST);
  switch (c) {
    default: m_progressiveAdd(&camera.par[c][0], d); break;
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
  FILE *f;

  int width = config.width;
  int height = config.height;

  if ((f = fopen(tgaFile, "wb")) != 0) {
    unsigned char header[18] = {
      0,0,2,0,0,0,0,0,0,0,0,0,width%256,width/256,height%256,height/256,24,0
    };
    unsigned char* img = (unsigned char*)malloc(width * height * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_FRONT);
    glReadPixels(viewportOffset[0], viewportOffset[1], width, height, GL_BGR, GL_UNSIGNED_BYTE, img);

    fwrite(header, 18, 1, f);
    fwrite(img, 3, width*height, f);

    free(img);
    fclose(f);
    printf(__FUNCTION__ " : wrote %s\n", tgaFile);
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

// Read out Z-buffer around the center of the view.
float distanceToSurface() {
#if TRY_DISTANCE_TO_SURFACE
  const int SIZE = 1;
  float z[SIZE*SIZE];
  int x = config.width / 2 - SIZE / 2;
  int y = config.height / 2 - SIZE / 2;
  glReadPixels(viewportOffset[0] + x, viewportOffset[1] + config.height - 1 - y, SIZE, SIZE,
               GL_DEPTH_COMPONENT, GL_FLOAT, z);
  float avg = 0;
  for (size_t i = 0; i < lengthof(z); ++i) {
    avg += z[i];
  }
  float v = avg / lengthof(z);

  // Convert zbuffer [0,1] value v back into actual dist.
  const float a = zFar / (zFar - zNear);
  const float b = zFar * zNear / (zNear - zFar);

  return b / (v - a);
#else
  return 100.0f;
#endif
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

  glsl_source.assign(fs);

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

  if (vs != default_vs) free((char*)vs);
  if (fs != default_fs) free((char*)fs);

  if (glGetError()) die("setupShaders() fails");

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

  if (glGetError()) die("setupShaders() fails");

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

  glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                   GL_LINEAR_MIPMAP_LINEAR);

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
    char xyz = 'x';
    if (sscanf(line.c_str() + parStart + 5, "%d].%c",
               &index, &xyz) < 1) continue;
    if (index < 0 || index > (int)lengthof(camera.par)) continue;
    if (xyz < 'x' || xyz > 'z') continue;

    printf("parameter %s par[%d].%c {%s}\n",
           varName.c_str(), index, xyz, attr.c_str());

    float* address = &camera.par[index][xyz-'x'];

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
    sprintf(filename, "keyframe-%u.cfg", i);
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
    sprintf(filename, "keyframe-%lu.cfg", (unsigned long)i);
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
  // Peel known options off the back..
  while (argc>1) {
    if (!strcmp(argv[argc-1], "--overunder")) {
      stereoMode = ST_OVERUNDER;
    } else if (!strcmp(argv[argc-1], "--interlaced")) {
      stereoMode = ST_INTERLACED;
    } else if (!strcmp(argv[argc-1], "--xeyed")) {
      stereoMode = ST_XEYED;
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
    } else break;
    --argc;
  }

  const char* configFile = (argc>=2 ? argv[1] : DEFAULT_CONFIG_FILE);

  // Load configuration.
  if (config.loadConfig(configFile)) {
    changeWorkingDirectory(configFile);
  } else {
    die("Usage: boxplorer <configuration-file.cfg>");
  }

  // Sanitize / override config parameters.
  if (loop) config.loop = true;
  if (enableDof) config.enable_dof = (enableDof == 1);  // override
  if (config.fps < 5) config.fps = 30;
  if (config.depth_size < 16) config.depth_size = 16;
  if (stereoMode == ST_XEYED) config.width *= 2; 

  camera = config;

  bool keyframesChanged = false;
  LoadKeyFrames(fixedFov);

  // Initialize SDL and OpenGL graphics.
  SDL_Init(SDL_INIT_VIDEO) == 0 ||
      die("SDL initialization failed: %s\n", SDL_GetError());
  atexit(SDL_Quit);

  // Set up the video mode, OpenGL state, shaders and shader parameters.
  initGraphics();
  initTwBar();
  initFPS(FPS_FRAMES_TO_AVERAGE);

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

    if (!rendering) {
      // If we're rendering a sequence to disk, we don't care about z-buffer.
      // Otherwise, just overwrite since we write every pixel.
      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_ALWAYS);
    }

    glUseProgram(program);
    camera.render(stereoMode);
    glUseProgram(0);

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
        glVertex3fv(splines[i].pos());
        splines[i].move(-2*splines[i].speed, 0, 0);
      }
      glEnd();
      glBegin(config.loop?GL_LINE_LOOP:GL_LINE_STRIP);  // left eye
      for (size_t i = 0; i < splines.size(); ++i) {
        glVertex3fv(splines[i].pos());
      }
      glEnd();

      glLineWidth(13);
      for (size_t i = 0; i < keyframes.size(); ++i) {
        glColor4f(0,0,1 - (i/256.0),1);  // Encode keyframe # in color.
        glBegin(GL_LINES);
        KeyFrame tmp = keyframes[i];
        tmp.move(tmp.speed, 0, 0);
        glVertex3fv(tmp.pos());
        tmp = keyframes[i];
        tmp.move(-tmp.speed, 0, 0);
        glVertex3fv(tmp.pos());
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

    float dist = distanceToSurface();

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
    sprintf(caption, "%s %.2ffps %5lu [%.3f %.3f %.3f] %dms",
        printController(controllerStr, ctl),
        getFPS(), (unsigned long)splines_index,
        camera.pos()[0], camera.pos()[1], camera.pos()[2],
        getLastFrameDuration()
    );
    SDL_WM_SetCaption(caption, 0);

    // Process events.
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (grabbedInput || !TwEventSDL(&event)) switch (event.type) {
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
          case 4:  // mouse wheel up, increase speed at keyframe
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].speed *= 1.1;
             }
            break;
          case 5:  // mouse wheel down, decrease speed at keyframe
            if (keyframe < keyframes.size()) {
              keyframes[keyframe].speed *= .9;
             }
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

      case SDL_KEYDOWN: switch (event.key.keysym.sym) {
       case SDLK_ESCAPE: {
         if (grabbedInput && !config.fullscreen) {
           grabbedInput = 0;
           SDL_ShowCursor(SDL_ENABLE);
           SDL_WM_GrabInput(SDL_GRAB_OFF);
         } else done |= 1;
      } break;

      // Switch fullscreen mode (loses the whole OpenGL context in Windows).
      case SDLK_RETURN: case SDLK_KP_ENTER: {
        config.fullscreen ^= 1; grabbedInput = 1; initGraphics(); initTwBar();
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
          if (event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL)) {
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
        if (event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL)) {
          // Replace currently selected keyframe.
          index = keyframe;
        }
        if (!(event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL))) {
          // Need an estimate for delta_time for this new frame.
          suggestDeltaTime(camera, keyframes);
          keyframes.push_back(camera);  // Add keyframe at end.
        } else {
          keyframes[index] = camera;  // Overwrite current keyframe.
        }
        char filename[256];
        sprintf(filename, "keyframe-%lu.cfg", (unsigned long)index);
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
        if (!(event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL))) ++keyframe;
        if (keyframe >= keyframes.size()) keyframe = 0;

        if (keyframe < keyframes.size() && splines.empty()) {
           // Don't jump camera ahead if we were playing, just stop in place.
           camera = keyframes[keyframe];
        }

        if (keyframe < keyframes.size()) {
          printf("at keyframe %lu, speed %f, delta_time %f\n",
              (unsigned long)keyframe, keyframes[keyframe].speed,
              keyframes[keyframe].delta_time);
        }

        if (keyframes.empty()) camera = config;  // back to start

        if (event.key.keysym.mod & (KMOD_LCTRL|KMOD_RCTRL)) {
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
          printf("at keyframe %lu, speed %f, delta_time %f\n",
              (unsigned long)keyframe, keyframes[keyframe].speed,
              keyframes[keyframe].delta_time);
        } else camera = config;
        splines.clear();
      } break;

      case SDLK_LSHIFT: {
        // Change movement speed.
        if (camera.speed < 1) camera.speed *= 1.1;
      } break;

      case SDLK_LALT: {
        if (camera.speed > camera.min_dist) camera.speed *= .9;
      } break;

      // Resolve controller value changes that happened during rendering.
      case SDLK_LEFT:  ctlXChanged = 1; updateControllerX(ctl, -(consecutiveChanges=1)); break;
      case SDLK_RIGHT: ctlXChanged = 1; updateControllerX(ctl,  (consecutiveChanges=1)); break;
      case SDLK_DOWN:  ctlYChanged = 1; updateControllerY(ctl, -(consecutiveChanges=1)); break;
      case SDLK_UP:    ctlYChanged = 1; updateControllerY(ctl,  (consecutiveChanges=1)); break;

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
      }
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

    // Continue after calling SDL_GetRelativeMouseState() so view direction
    // does not jump after closing AntTweakBar.
    if (!grabbedInput) continue;

    (void)mouse_buttons;
    (void)mouse_button_left;
    (void)mouse_button_right;

    if (keystate[SDLK_w]) camera.move(0, 0, min(camera.speed , 0.9f * dist));  //forward
    if (keystate[SDLK_s]) camera.move(0, 0, -camera.speed);  //back

    if (keystate[SDLK_a]) camera.move(-camera.speed, 0, 0);  //left
    if (keystate[SDLK_d]) camera.move( camera.speed, 0, 0);  //right

    // Mouse look.
    if (grabbedInput && (mouse_dx != 0 || mouse_dy != 0)) {
      m_rotateX2(camera.mouse_rot_speed * mouse_dx);
      m_rotateY2(camera.mouse_rot_speed * mouse_dy);
    }
    if (keystate[SDLK_q]) m_rotateZ2(camera.keyb_rot_speed);
    if (keystate[SDLK_e]) m_rotateZ2(-camera.keyb_rot_speed);

   if (keystate[SDLK_z]){ if (camera.speed > camera.min_dist) camera.speed -= camera.speed/10; printf("speed %f\n", camera.speed);}
   if (keystate[SDLK_c]){ if (camera.speed < 1.0) camera.speed += camera.speed/10; printf("speed %f\n", camera.speed);}

    // Change the value of the active controller.
    if (!ctlXChanged) {
      if (keystate[SDLK_LEFT])  { ctlXChanged = 1; updateControllerX(ctl, -++consecutiveChanges); }
      if (keystate[SDLK_RIGHT]) { ctlXChanged = 1; updateControllerX(ctl,  ++consecutiveChanges); }
    }
    if (!ctlYChanged) {
      if (keystate[SDLK_DOWN])  { ctlYChanged = 1; updateControllerY(ctl, -++consecutiveChanges); }
      if (keystate[SDLK_UP])    { ctlYChanged = 1; updateControllerY(ctl,  ++consecutiveChanges); }
    }

    if (!(ctlXChanged || ctlYChanged)) consecutiveChanges = 0;
  }

  TwTerminate();

  if (!rendering && keyframesChanged) {
    // TODO: ask whether to save keyframes
  }
  camera.saveConfig("last.cfg");  // Save a config file on exit, just in case.
  return 0;
}
