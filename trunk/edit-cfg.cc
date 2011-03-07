// Little cmdline utility to change or delete values from .cfg files.
// Useful for batch editing of set of keyframes.
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;

int main(int argc, char* argv[]) {
  int arg = 1;
  bool doDelete = false;

  if (arg < argc && !strcmp(argv[arg], "-d")) {
    doDelete = true;
    ++arg;
  }
  string kv(arg<argc?argv[arg++]:"");

  if (kv.empty()) {
    fprintf(stderr, "usage: edit-cfg [-d] \"param value\" *.cfg\n");
    return 1;
  }

  string key(kv);
  if (key.find(' ') != string::npos) key.erase(key.find(' ') + 1);

  // arg points to file(s) to process
  while (arg < argc) {
    string filename(argv[arg++]);
    printf("processing %s: ", filename.c_str());

    vector<string> content;

    // read file content
    ifstream infile(filename.c_str());
    if (!infile.is_open()) {
      printf("failed to open for read\n");
      continue;
    }
    string line;
    while(infile.good()) {
      getline(infile,line);
      content.push_back(line);
    }
    infile.close();

    // change key value
    for (vector<string>::iterator it = content.begin();
         it != content.end(); ++it) {
      if (it->find(key) == 0) {
        it = content.erase(it);
        if (!doDelete) it = content.insert(it, kv);
      }
    }

    // write file
    ofstream outfile(filename.c_str());
    if (outfile.is_open()) {
      for(vector<string>::const_iterator it = content.begin();
          it != content.end(); ++it) {
        outfile.write(it->data(), it->size());
        outfile.write("\r\n", 2);
      }
      outfile.close();
    }
    printf("done\n");
  }

  return 0;
}
