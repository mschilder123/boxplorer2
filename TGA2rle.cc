#include "TGA.h"
#include <stdio.h>
#include <stdlib.h>

const int patterns[] = {
    0x3030, 0x6060, 0xc0c0, 0x0303, 0x0606, 0x0c0c,

    0x0099, 0x0990, 0x9900,

    0x0055, 0x0550, 0x5500,

    0x0306, 0x0603, 0x3060, 0x6030,

    0x060c, 0x0c06, 0x60c0, 0xc060,

    0x2980,

    0x0154, 0x1540, 0x02a8, 0x2a80,
};

int main(int argc, char *argv[]) {
  TGA tga;
  tga.readFile(argv[1]);

  for (int y = 0; y < tga.height(); y += 4) {
    for (int x = 0; x < tga.width(); x += 4) {
      bool set = false;
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          set = set || tga.get(x + i, y + j);
          tga.set(x + i, y + j, 0, 0, 0);
        }
      }
      if (set) {
        // pick random 4x4 walker pattern
        int pattern =
            patterns[rand() % (sizeof(patterns) / sizeof(patterns[0]))];
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            if (pattern & (1 << (i + 4 * j))) {
              tga.set(x + i, y + j, 255, 255, 255);
            }
          }
        }
      }
    }
  }

  printf("x = %d, y = %d, rule = SingleRotation\n", tga.width(), tga.height());

  int prev_y = 0;
  bool empty_line = true;
  for (int y = 0; y < tga.height(); ++y) {
    bool prev = tga.get(0, tga.height() - 1 - y);
    int prev_x = 0;
    for (int x = 1; x < tga.width(); ++x) {
      bool cur = tga.get(x, tga.height() - 1 - y);
      if (cur != prev) {
        if (x - prev_x > 1)
          printf("%d", x - prev_x);
        printf("%c", prev ? 'o' : 'b');
        prev = cur;
        prev_x = x;
        empty_line = false;
      }
    }
    if (empty_line == true && prev) { // full line
      printf("%d%c", tga.width(), 'o');
      empty_line = false;
    }
    if (empty_line == false || y == tga.height() - 1) {
      if (y - prev_y > 1)
        printf("%d", y - prev_y);
      printf("$\n");
      prev_y = y;
      empty_line = true;
    }
  }
  printf("!");

  return 0;
}
