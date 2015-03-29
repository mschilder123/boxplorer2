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
  bool doRename = false;
  bool doClean = false;

  if (arg < argc && argv[arg][0] == '-') {
    if (!strcmp(argv[arg], "-d"))
      doDelete = true;
    if (!strcmp(argv[arg], "-r"))
      doRename = true;
    if (!strcmp(argv[arg], "-c"))
      doClean = true;
    ++arg;
  }

  string kv(arg<argc?argv[arg++]:"");

  if (!doClean && kv.empty()) {
    fprintf(stderr, "usage: edit-cfg [-drc] \"param value\" *.cfg\n");
    return 1;
  }

  string key(kv);
  string value;

  size_t space = key.find(' ');
  if (space != string::npos) {
    value.assign(key.substr(space + 1));
    key.erase(space + 1);
  }

  fprintf(stderr, "key %s, value %s\n", key.c_str(), value.c_str());

  // arg points to file(s) to process
  bool cleaning = false;
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
      getline(infile, line);
      if (line.find_first_of("\r\n") != string::npos)
        line.erase(line.find_first_of("\r\n"));
      // Clean svn diff pollution..
      if (doClean && line.find("<<<<<") == 0) {
              cleaning = true;
              continue;
      }
      if (doClean && line.find("=====") == 0) {
              cleaning = false;
              continue;
      }
      if (doClean && line.find(">>>>>") == 0) {
              cleaning = false;
              continue;
      }
      if (!cleaning && !line.empty())
        content.push_back(line);
    }
    infile.close();

    if (!doClean) {
      // change key value
      bool found = false;

      for (vector<string>::iterator it = content.begin();
           it != content.end(); ++it) {
        if (doRename && it->find(value) == 0) {
          string line(*it);
          it->clear();
        } else
        if (it->find(key) == 0) {
          string line(*it);
          it->clear();
          if (!doDelete && !doRename) it->assign(kv);
          if (doRename) {
            string newline = line.replace(0, key.size() - 1, value);
            it->assign(newline);
          }
          found = true;
        }
      }

      if (!found && !doDelete && !doRename)
        content.insert(content.begin(), kv);
    }

    // write file
    ofstream outfile(filename.c_str());
    if (outfile.is_open()) {
      for(vector<string>::const_iterator it = content.begin();
          it != content.end(); ++it) {
        if (!it->empty()) {
          outfile.write(it->data(), it->size());
          outfile.write("\n", 1);
        }
      }
      outfile.close();
    }
    printf("done\n");
  }

  return 0;
}
