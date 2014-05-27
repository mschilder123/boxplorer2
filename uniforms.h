#ifndef _F_UNIFORMS_H_
#define _F_UNIFORMS_H_

#include <string>

#if defined(__APPLE__) || defined(__GNUC__)
#include <unordered_map>
#define hash_map unordered_map
#else  // WIN32
#include <hash_map>
#endif

class KeyFrame;

// Generic uniform interface.
class iUniform {
 public:
  virtual ~iUniform(){}
  virtual const std::string& name() = 0;
  virtual std::string toString() = 0;
  virtual bool parse(const std::string& line) = 0;
  virtual void send(int prog) = 0;  // send to shader
  virtual void bindToUI(void* bar) = 0;  // register as tw var
  virtual bool link(KeyFrame* kf) = 0;  // link to storage within keyframe.
  virtual iUniform* Clone() = 0;
  virtual bool ok() = 0;
};

// Uniform interface instance smart ptr.
// Can be used safely in hash_map<>
class iUniformPtr {
 public:
  iUniformPtr() : ptr_(NULL) {}
  explicit iUniformPtr(iUniform* ptr) : ptr_(ptr) {}
  // Copy constructor needs to Clone() ptr_.
  iUniformPtr(const iUniformPtr& other) : ptr_(NULL) {
    if (other.ptr_) ptr_ = other.ptr_->Clone();
  }
  // Assignment operator also needs to Clone() ptr_.
  iUniformPtr& operator=(const iUniformPtr& other) {
    if (ptr_) { delete ptr_; ptr_ = NULL; }
    if (other.ptr_) ptr_ = other.ptr_->Clone();
    return *this;
  }
  iUniform* operator->(void) { return ptr_; }
  ~iUniformPtr() { delete ptr_; }
  bool ok() { return ptr_ != NULL && ptr_->ok(); }
 private:
  bool operator==(const iUniformPtr& other);
  bool operator<(const iUniformPtr& other);
  iUniform* ptr_;
};

// Our active uniforms.
class Uniforms {
 public:
  // Try parse declared uniforms from GLSL source,
  // including UI attributes from comment.
  bool parseFromGlsl(const std::string& glsl);

  // Map a uniform var to local storage backed by kf.
  void link(KeyFrame* kf);

  // Register linked vars with UI.
  void bindToUI(void* bar);

  // Send to shader.
  void send(int program);

private:
  bool parseLine(const std::string& line, iUniformPtr* uni);

#if defined(__GNUC__) || defined(__APPLE__)
  std::hash_map<std::string, iUniformPtr, std::hash<std::string> > uniforms;
#else
  std::hash_map<std::string, iUniformPtr> uniforms;
#endif
};

#endif  // F_UNIFORMS_H_

