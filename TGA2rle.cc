#include "TGA.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  TGA tga;
  tga.readFile(argv[1]);

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
