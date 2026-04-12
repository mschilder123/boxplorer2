#include "gfx.h"

#include "shader_procs.h"
#include "utils.h"

// Defined in boxplorer2.cc; releases GL resources before context destruction.
void clearGlContext();

GFX::GFX()
    : display_(-1), width_(0), height_(0), last_x_(SDL_WINDOWPOS_CENTERED),
      last_y_(SDL_WINDOWPOS_CENTERED), last_width_(0), last_height_(0),
      fullscreen_(false), window_(NULL), glcontext_(0) {}

GFX::~GFX() { reset(); }

void GFX::init() {
  // requires SDL_Init(SDL_INIT_VIDEO);
  for (ndisplays_ = 0;
       ndisplays_ < SDL_GetNumVideoDisplays() && ndisplays_ < kMAXDISPLAYS;
       ++ndisplays_) {
    SDL_GetCurrentDisplayMode(ndisplays_, &mode_[ndisplays_]);
    SDL_GetDisplayBounds(ndisplays_, &rect_[ndisplays_]);
  }
}

void GFX::reset() {
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

bool GFX::resize(int w, int h) {
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

  DEBUG("%dx%d display %d", w, h, d);
  SDL_GL_SetAttribute(SDL_GL_FRAMEBUFFER_SRGB_CAPABLE, 1);
  window_ = SDL_CreateWindow("test", last_x_, last_y_, w, h,
                             SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
  glcontext_ = SDL_GL_CreateContext(window_);
  display_ = d;
  last_width_ = width_ = w;
  last_height_ = height_ = h;
  return enableShaderProcs(); // re-fetch ptrs in new context.
}

bool GFX::toggleFullscreen() {
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
    DEBUG("screen %d: %dx%d", i, rect_[i].w, rect_[i].h);
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
    DEBUG("to fullscreen %dx%d display %d", targetWidth, targetHeight, d);
    window_ = SDL_CreateWindow(
        "boxplorer2", rect_[d].x, rect_[d].y, targetWidth, targetHeight,
        SDL_WINDOW_OPENGL | SDL_WINDOW_FULLSCREEN_DESKTOP);
    width_ = targetWidth;
    height_ = targetHeight;
  } else {
    DEBUG("from fullscreen %dx%d display %d", last_width_, last_height_, d);
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
