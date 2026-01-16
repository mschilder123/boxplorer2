#include "utils.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#if !defined(_WIN32)
#include <unistd.h>
#else
#pragma warning(disable : 4996) // unsafe function
#endif

using namespace std;

string BaseDir;
string WorkingDir;
string BaseFile;

// Allocate a char[] and read a text file into it. Return 0 on error.
char *_readFile(char const *name) {
  FILE *f;
  size_t len;
  char *s = 0;

  // open file an get its length
  if (!(f = fopen(name, "r")))
    goto readFileError1;
  fseek(f, 0, SEEK_END);
  len = ftell(f);

  // read the file in an allocated buffer
  if (!(s = (char *)malloc(len + 1)))
    goto readFileError2;
  rewind(f);
  len = fread(s, 1, len, f);
  s[len] = '\0';

readFileError2:
  fclose(f);
readFileError1:
  return s;
}

bool readFile(const string &name, string *content) {
  string filename(WorkingDir + name);
  char *s = _readFile(filename.c_str());
  if (!s) {
    filename = BaseDir + "include/" + name;
    s = _readFile(filename.c_str());
    if (!s)
      return false;
  }
  content->assign(s);
  free(s);
  printf("%s : read '%s'\n", __func__, filename.c_str());
  return true;
}
