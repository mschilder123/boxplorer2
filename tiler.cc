// Simple reflective tiler, taking input image and
// mirroring it four ways into a tile.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#if !defined(__FUNCTION__)
#define __FUNCTION__ "tiler"
#endif

#include "TGA.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 2) {
     cerr << "Usage: " << argv[0] << " <tga-infile> [<tga-outfile>]" << endl;
     exit(1);
  }

  TGA input;
  if (!input.readFile(argv[1])) {
     cerr << "Failed to read '" << argv[1] << "'" << endl;
     exit(1);
  }

  TGA output(input.width() * 2, input.height() * 2);

  for (size_t x = 0; x < input.width(); ++x) {
    for (size_t y = 0; y < input.width(); ++y) {
      for (size_t c = 0; c < 3; ++c) {
        unsigned v = input.data()[x * 3 + c + y * input.width() * 3];
        // 0,0
        output.data()[x * 3 + c + y * output.width() * 3] = v;
        // 1,0
        output.data()[input.width() * 2 * 3 - 3 - x * 3 + c + y * output.width() * 3] = v;
        // 0,1
        output.data()[x * 3 + c + ((-y - 1) * output.width() + output.width() * input.height() * 2) * 3] = v;
        // 1,1
        output.data()[input.width() * 2 * 3 - 3 - x * 3 + c + ((-y - 1)* output.width() + output.width() * input.height() * 2) * 3] = v;
      }
    }
  }

  const char* outName = (argc > 2 ? argv[2] : argv[1]);
  if (!output.writeFile(outName)) {
     cerr << "Failed to write '" << outName << "'" << endl;
     exit(1);
  }

  return 0;
}
