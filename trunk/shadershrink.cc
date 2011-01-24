// Shrink a shader file and put it in a C header as a string.

#include <stdio.h>
#include <string.h>
#include <ctype.h>

int in_comment = 0;
char last_char = '\n';
int brace_depth = 1;
int anything_written = 0;

// Strip a shader line from spaces and comments and print it.
void print_stripped_line(char* s) {
  int preprocessor_line = 0;
  int i;
  int current_brace_depth = brace_depth;
  int line_start = 1;

  for (i=0; s[i]!='\n' && s[i]!=0; i++) {
    if (in_comment) {
      // end of "/* */" comment? treat it like a space
      if (s[i]=='*' && s[i+1]=='/') { in_comment = 0; s[i+1]=' '; continue; }
    }
    else {
      // treat "/* */" like spaces
      if (s[i]=='/' && s[i+1]=='*') { in_comment = 1; i++; continue; }
      // remove "//" comments
      else if (s[i]=='/' && s[i+1]=='/') { break; }
      // remove spaces that are not between alphanumeric characters
      else if (s[i]==' ' && !(isalnum(last_char) && isalnum(s[i+1]))) { continue; }
      // put a newline after preprocessor directives
      else if (s[i]=='#') { preprocessor_line = 1; }

      if (line_start) {
        line_start = 0;
        anything_written = 1;

        int spaces = brace_depth*2;
        if (last_char!='{' && last_char!='}' && last_char!='\n' && last_char!=';') spaces += 2;
        if (s[i]=='{' || s[i]=='}') spaces -= 2;
        if (s[i]==' ') spaces--;

        putchar('\n');
        while (spaces--) putchar(' ');
        putchar('"');
      }

      last_char = s[i];
      putchar(s[i]);
      if (s[i]=='"') putchar('"');

      if (s[i]=='{') current_brace_depth++;
      if (s[i]=='}') current_brace_depth--;
    }
  }

  if (preprocessor_line) { printf("\\n"); last_char = '\n'; }
  else { brace_depth = current_brace_depth; }

  if (!line_start) putchar('"');
}


int main(int argc, char** argv) {
  char s[32768] = {' '};  // insert a space at the beginning of every line

  if (argc < 2) return -1;
  printf("const char %s[] = ", argv[1]);

  while (fgets(s+1, sizeof(s)-1, stdin)) print_stripped_line(s);
  if (!anything_written) printf("\"\"");
  printf(";\n\n");

  return 0;
}
