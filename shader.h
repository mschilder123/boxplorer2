#ifndef _F_SHADER_H_
#define _F_SHADER_H_

#include <string>
#include <unordered_map>

class Shader {
public:
  Shader() : program_(0), ok_(false) {}
  virtual ~Shader() { clear(); }

  // Create, load and link shader.
  bool compile(const std::string &defines, const std::string &vertex_shader,
               const std::string &pixel_shader);

  // Detach and delete shader.
  void clear();

  unsigned int program() { return program_; }
  bool ok() const { return ok_; }

  const std::string &source() const { return source_; }
  const std::string &log() const { return log_; }
  const std::string &uniforms() const { return uniforms_; }

  int uniform_location(const std::string &name) const {
    auto it = uniform_locations_.find(name);
    if (it == uniform_locations_.end())
      return -1;
    return it->second;
  }

private:
  std::string source_;
  unsigned int program_;
  std::string log_;
  std::string uniforms_;
  std::unordered_map<std::string, int> uniform_locations_;
  bool ok_;
};

#endif // _F_SHADER_H_
