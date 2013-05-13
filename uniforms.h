#ifndef _F_UNIFORMS_H_
#define _F_UNIFORMS_H_

#include <string>

class KeyFrame;

// Generic uniform interface.
class iUniform {
 public:
  virtual ~iUniform(){}
  virtual const std::string& name() = 0;
  virtual std::string toString() = 0;
  virtual void send(int prog) = 0;  // send to shader
  virtual void twVar(void* bar) = 0;  // register tw var
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

// Map a uniform var to local storage backed by kf.
iUniformPtr link_uniform(const std::string& line, KeyFrame* kf);

#endif
