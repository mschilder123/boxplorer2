
#pragma once
#include <string>

// Global variables for paths
extern std::string BaseDir;
extern std::string WorkingDir;
extern std::string BaseFile;

// Utility functions
bool readFile(const std::string &name, std::string *content);
char *_readFile(char const *name); // Internal helper exposed if needed
