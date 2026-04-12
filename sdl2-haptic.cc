#include <SDL.h>
#include <iostream>

#if defined(_WIN32)
#pragma comment(lib, "SDL2.lib")
#pragma comment(lib, "SDL2main.lib")
#pragma comment(lib, "Shell32.lib")
// SDL2 statics
#pragma comment(lib, "advapi32.lib")
#pragma comment(lib, "version.lib")
#pragma comment(lib, "imm32.lib")
#pragma comment(lib, "setupapi.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "ole32.lib")
#pragma comment(lib, "oleaut32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "winmm.lib")
#endif

int main(int argc, char *argv[]) {
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK | SDL_INIT_HAPTIC) < 0) {
    std::cerr << "SDL could not initialize! SDL Error: " << SDL_GetError()
              << std::endl;
    return 1;
  }

  SDL_Joystick *joystick = nullptr;
  SDL_Haptic *haptic = nullptr;

  // 1. Check for joysticks
  if (SDL_NumJoysticks() < 1) {
    std::cerr << "No joysticks connected!" << std::endl;
  } else {
    // 2. Open the first joystick
    joystick = SDL_JoystickOpen(0);
    if (joystick == nullptr) {
      std::cerr << "Warning: Unable to open joystick! SDL Error: "
                << SDL_GetError() << std::endl;
    } else {
      // 3. Check if it has haptic capability
      if (SDL_JoystickIsHaptic(joystick)) {
        haptic = SDL_HapticOpenFromJoystick(joystick);
        if (haptic == nullptr) {
          std::cerr << "Warning: Unable to open haptic! SDL Error: "
                    << SDL_GetError() << std::endl;
        }
      }
    }
  }

  if (haptic != nullptr) {
    std::cout << "Haptic device opened: " << SDL_HapticName(0) << std::endl;

    // 4. Initialize simple rumble or complex effects
    if (SDL_HapticRumbleSupported(haptic)) {
      SDL_HapticRumbleInit(haptic);
    }

    // --- Custom FFB Effect (Constant Force) ---
    SDL_HapticEffect effect;
    SDL_memset(&effect, 0, sizeof(SDL_HapticEffect));
    // effect.type = SDL_HAPTIC_CONSTANT;
    effect.type = SDL_HAPTIC_FRICTION;
    effect.constant.direction.type = SDL_HAPTIC_POLAR; // Polar coordinates
    effect.constant.direction.dir[0] = 0;              // 0 degrees = right
    effect.constant.length = 5000;        // Duration in ms (5 seconds)
    effect.constant.level = 20000;        // Strength (-32767 to 32767)
    effect.constant.attack_length = 1000; // Fade in time

    // 5. Upload the effect
    int effect_id = SDL_HapticNewEffect(haptic, &effect);

    // 6. Run the effect
    SDL_HapticRunEffect(haptic, effect_id, 1);
    std::cout << "Running constant force effect..." << effect_id << std::endl;

    SDL_Delay(5000); // Wait 5 seconds

    // 7. Clean up
    SDL_HapticDestroyEffect(haptic, effect_id);
  }

  // Cleanup
  if (haptic != nullptr)
    SDL_HapticClose(haptic);
  if (joystick != nullptr)
    SDL_JoystickClose(joystick);
  SDL_Quit();
  return 0;
}
