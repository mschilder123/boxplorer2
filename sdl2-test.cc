// SDL2 sample showing how to resize and toggle fullscreen.
//
// You can resize and drag the window across screens and
// <ENTER> will toggle fullscreen on its current display.
//
// Handy for multi-monitor setups like passive 3d secondaries or oculus.
//
// marius.schilder@gmail.com

#include <stdio.h>

#include <SDL.h>
#include <SDL_main.h>

#if defined(_WIN32)
#pragma comment(lib, "SDL2.lib")
#pragma comment(lib, "SDL2main.lib")
#endif

#if !defined(__FUNCTION__)
#define __FUNCTION__ "sdl2-test"
#endif

#define MAXDISPLAYS 6

class GFX {
 public:
  GFX() : display_(-1),
          width_(0), height_(0),
          last_x_(SDL_WINDOWPOS_CENTERED),
          last_y_(SDL_WINDOWPOS_CENTERED),
          last_width_(0), last_height_(0),
          fullscreen_(false),
          window_(NULL), renderer_(NULL) {
    SDL_Init(SDL_INIT_VIDEO);
    for (int i = 0; i < SDL_GetNumVideoDisplays() &&
                    i < MAXDISPLAYS; ++i) {
      SDL_GetCurrentDisplayMode(i, &mode_[i]);
      SDL_GetDisplayBounds(i, &rect_[i]);
    }
  }

  ~GFX() {
    reset();
    SDL_Quit();
  }

  void reset() {
    SDL_DestroyRenderer(renderer_);
    renderer_ = NULL;
    SDL_DestroyWindow(window_);
    window_ = NULL;
    display_ = -1;
  }

  void resize(int w, int h) {
    int d = display_;

    // ignore resize events when fullscreen.
    if (fullscreen_) return;

    // ignore resize events for fullscreen width.
    if (d != -1 && w == rect_[d].w) return;

    if (window_) {
      // capture current display.
      d = SDL_GetWindowDisplayIndex(window_);
      // capture current window position.
      SDL_GetWindowPosition(window_, &last_x_, &last_y_);
    }

    reset();

    printf(__FUNCTION__ ": %dx%d display %d\n", w, h, d);
    window_ = SDL_CreateWindow("test",
       last_x_, last_y_,
       w, h,
       SDL_WINDOW_OPENGL|SDL_WINDOW_BORDERLESS);
    renderer_ = SDL_CreateRenderer(window_, d, 0);
    display_ = d;
    last_width_ = width_ = w;
    last_height_ = height_ = h;
  }

  void toggleFullscreen() {
    if (!window_) return;

    // capture current display.
    int d = SDL_GetWindowDisplayIndex(window_);

    if (!fullscreen_) {
      // capture current window position.
      SDL_GetWindowPosition(window_, &last_x_, &last_y_);
    }
    reset();

    if (!fullscreen_) {
      printf(__FUNCTION__ ": to fullscreen %dx%d display %d\n",
                      rect_[d].w, rect_[d].h, d);
      window_ = SDL_CreateWindow("test",
          rect_[d].x, rect_[d].y,
          rect_[d].w, rect_[d].h,
          SDL_WINDOW_OPENGL|SDL_WINDOW_FULLSCREEN);
      width_ = rect_[d].w;
      height_ = rect_[d].h;
    } else {
      printf(__FUNCTION__ ": from fullscreen %dx%d display %d\n",
                      last_width_, last_height_, d);
      window_ = SDL_CreateWindow("test",
          last_x_, last_y_,
          last_width_, last_height_,
          SDL_WINDOW_OPENGL|SDL_WINDOW_BORDERLESS);
      width_ = last_width_;
      height_ = last_height_;
    }

    renderer_ = SDL_CreateRenderer(window_, d, 0);
    display_ = d;

    fullscreen_ = !fullscreen_;
  }

  SDL_Renderer* renderer() { return renderer_; }
  int width() const { return width_; }
  int height() const { return height_; }

 private:
  int display_;
  int width_, height_;  // current dimensions, window or fullscreen.
  int last_x_,last_y_;  // last known position of window
  int last_width_, last_height_;  // last known dimension of window.
  bool fullscreen_;
  SDL_Window* window_;
  SDL_Renderer* renderer_;
  SDL_DisplayMode mode_[MAXDISPLAYS];
  SDL_Rect rect_[MAXDISPLAYS];
};

int main(int argc, char* argv[]) {
  GFX gfx;

  gfx.resize(720, 480);

  int frame = 0;
  bool done = false;

  SDL_InitSubSystem(SDL_INIT_JOYSTICK);
  SDL_Joystick* joystick = SDL_JoystickOpen(2);

    printf(__FUNCTION__ " : JoystickName '%s'\n",
             SDL_JoystickName(joystick));
    printf(__FUNCTION__ " : JoystickNumAxes   : %i\n",
           SDL_JoystickNumAxes(joystick));
    printf(__FUNCTION__ " : JoystickNumButtons: %i\n",
           SDL_JoystickNumButtons(joystick));
    printf(__FUNCTION__ " : JoystickNumHats   : %i\n",
           SDL_JoystickNumHats(joystick));
   SDL_JoystickEventState(SDL_ENABLE);

  Sint16 axes[6] = {0};

  while(!done) {
    ++frame;

    // Draw some flickering background.
    SDL_SetRenderDrawColor(gfx.renderer(),
                    200*(frame&1), 100*(frame&2), 60*(frame*4), 255);

    // Clear the entire screen to our selected color.
    SDL_RenderClear(gfx.renderer());

    // Show updated render.
    SDL_RenderPresent(gfx.renderer());

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_KEYDOWN: {
          switch (event.key.keysym.sym) {
            case SDLK_RETURN: {
              gfx.toggleFullscreen();
            } break;
            case SDLK_ESCAPE: {
              done = true;
            } break;
          }
        } break;
        case SDL_WINDOWEVENT: {
          if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
            gfx.resize(event.window.data1, event.window.data2);
          }
        } break;
        case SDL_JOYBUTTONDOWN:
        case SDL_JOYBUTTONUP:
        {
          printf("%s:%d\n",
                 event.jbutton.type == SDL_JOYBUTTONDOWN?"dn":"up",
                 event.jbutton.button
          );
        } break;
        case SDL_JOYAXISMOTION:
        {
          Sint16 v = event.jaxis.value;
          if (v < -5000 || v > 5000) {
            axes[event.jaxis.axis] = v;
            for (int i = 0; i < 6; ++i) {
              printf("%d:%8d ",
                     i, axes[i]);
            }
            printf("\n");
          }
        } break;
      }
    }

    const Uint8* state =SDL_GetKeyboardState(NULL);
    if (state[SDL_SCANCODE_W]) {
            printf("W");
    }

    fflush(stdout);
  }

  return 0;
}
