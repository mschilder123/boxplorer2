#ifndef GFX_H_
#define GFX_H_

#define NO_SDL_GLEXT
#include <SDL.h>
#include <SDL_opengl.h>

static const int kMAXDISPLAYS = 6;

// Our OpenGL SDL2 window.
class GFX {
public:
  GFX();
  ~GFX();

  void init();
  void reset();
  bool resize(int w, int h);
  bool toggleFullscreen();

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
};

#endif // GFX_H_
