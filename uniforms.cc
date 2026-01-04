#include "uniforms.h"

#include <iostream>
#include <string>
#include <sstream>

#include <AntTweakBar.h>

#include "camera.h"
#include "shader_procs.h"

using namespace std;

namespace {
  // Parse a line into type,name,attr.
  // TODO: handle arrays?
  bool parseLine(const string& line,
               string* type,
               string* name,
               string* attr) {
    istringstream is(line);
    string uniform;
    is >> uniform;  // drop "uniform"
    is >> *type;
    is >> *name;
    size_t semi = name->find(';');
    if (semi == string::npos) return false;
    name->erase(semi);
    size_t at = line.find('{');
    if (at != string::npos) {
      size_t end_at = line.find('}', at);
      if (end_at != string::npos) {
        attr->assign(line, at + 1, end_at - at - 1);
      }
    }
    return true;
  }
}

class IntUniform : public iUniform {
 public:
   IntUniform() : adr_(NULL) {}
   virtual ~IntUniform() {}
   iUniform* Clone() { return new IntUniform(*this); }

   bool parse(const string& line) {
     return parseLine(line, &type_, &name_, &attr_);
   }
   const string& name() { return name_; }
   string toString() {
    ostringstream o;
    o << "int " << name_ << " ";
    if (adr_ != NULL) o << *adr_  << " @" << adr_;
    else o << "(nil)";
    return o.str();
   }
   void bindToUI(void* vbar) {
     TwAddVarRW((TwBar*)vbar, name_.c_str(), TW_TYPE_INT32, adr_,
                attr_.c_str());
   }
   void send(int program) {
     glUniform1i(glGetUniformLocation(program, name_.c_str()), *adr_);
   }
   bool link(KeyFrame* kf) {
     adr_ = (int*)kf->map_address(type_, name_, 1);
     return adr_ != NULL;
   }
   bool ok() { return adr_ != NULL; }
 private:
  IntUniform(const IntUniform& other) :
          adr_(other.adr_), name_(other.name_),
          type_(other.type_), attr_(other.attr_) {}
  IntUniform& operator=(const IntUniform& other);

  int* adr_;
  string name_;
  string type_;
  string attr_;
};

class BoolUniform : public iUniform {
 public:
   BoolUniform() : adr_(NULL) {}
   virtual ~BoolUniform() {}
   iUniform* Clone() { return new BoolUniform(*this); }

   bool parse(const string& line) {
     return parseLine(line, &type_, &name_, &attr_);
   }
   const string& name() { return name_; }
   string toString() {
    ostringstream o;
    o << "bool " << name_ << " ";
    if (adr_ != NULL) o << *adr_  << " @" << adr_;
    else o << "(nil)";
    return o.str();
   }
   void bindToUI(void* vbar) {
     TwAddVarRW((TwBar*)vbar, name_.c_str(), TW_TYPE_BOOL32, adr_,
                attr_.c_str());
   }
   void send(int program) {
     glUniform1i(glGetUniformLocation(program, name_.c_str()), *adr_);
   }
   bool link(KeyFrame* kf) {
     adr_ = (int*)kf->map_address("int", name_, 1);
     return adr_ != NULL;
   }
   bool ok() { return adr_ != NULL; }
 private:
  BoolUniform(const BoolUniform& other) :
          adr_(other.adr_), name_(other.name_),
          type_(other.type_), attr_(other.attr_) {}
  BoolUniform& operator=(const BoolUniform& other);

  int* adr_;
  string name_;
  string type_;
  string attr_;
};

class FloatUniform : public iUniform {
 public:
   FloatUniform() : adr_(NULL) {}
   virtual ~FloatUniform() {}
   iUniform* Clone() { return new FloatUniform(*this); }

   bool parse(const string& line) {
     return parseLine(line, &type_, &name_, &attr_);
   }
   const string& name() { return name_; }
   string toString() {
    ostringstream o;
    o << "float " << name_ << " ";
    if (adr_ != NULL) o << *adr_ << " @" << adr_;
    else o << "(nil)";
    return o.str();
   }
   void bindToUI(void* bar) {
     TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_FLOAT, adr_, attr_.c_str());
   }
   void send(int program) {
     glUniform1f(glGetUniformLocation(program, name_.c_str()), *adr_);
   }
   bool link(KeyFrame* kf) {
    adr_ = (float*)kf->map_address(type_, name_, 1);
    return adr_ != NULL;
   }
   bool ok() { return adr_ != NULL; }
private:
    FloatUniform(const FloatUniform& other) :
            adr_(other.adr_), name_(other.name_),
            type_(other.type_), attr_(other.attr_) {}
  FloatUniform& operator=(const FloatUniform& other);

  float* adr_;
  string name_;
  string type_;
  string attr_;
};

#if defined(GL_ARB_gpu_shader_fp64)
class DoubleUniform : public iUniform {
 public:
   DoubleUniform() : adr_(NULL) {}
   virtual ~DoubleUniform() {}
   iUniform* Clone() { return new DoubleUniform(*this); }

   bool parse(const string& line) {
     return parseLine(line, &type_, &name_, &attr_);
   }
   const string& name() { return name_; }
   string toString() {
    ostringstream o;
    o << "double " << name_ << " ";
    if (adr_ != NULL) o << *adr_ << " @" << adr_;
    else o << "(nil)";
    return o.str();
   }
   void bindToUI(void* bar) {
     TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_DOUBLE, adr_,
                attr_.c_str());
   }
   void send(int program) {
     glUniform1d(glGetUniformLocation(program, name_.c_str()), *adr_);
   }
   bool link(KeyFrame* kf) {
    adr_ = (double*)kf->map_address(type_, name_, 1);
    return adr_ != NULL;
   }
   bool ok() { return adr_ != NULL; }
private:
    DoubleUniform(const DoubleUniform& other) :
            adr_(other.adr_), name_(other.name_),
            type_(other.type_), attr_(other.attr_) {}
  DoubleUniform& operator=(const DoubleUniform& other);

  double* adr_;
  string name_;
  string type_;
  string attr_;
};
#endif

class Vec3Uniform : public iUniform {
 public:
   Vec3Uniform() : adr_(NULL) {}
   virtual ~Vec3Uniform() {}
   iUniform* Clone() { return new Vec3Uniform(*this); }

   bool parse(const string& line) {
    return parseLine(line, &type_, &name_, &attr_);
   }
   const string& name() { return name_; }
   string toString() {
    ostringstream o;
    o << "vec3 " << name_ << " ";
    if (adr_ != NULL) o << adr_[0] << " "
                        << adr_[1] << " " << adr_[2] << " @" << adr_;
    else o << "(nil)";
    return o.str();
   }
   void bindToUI(void* bar) {
     if (name_.find("Color") != string::npos) {
       TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_COLOR3F, adr_,
                       attr_.c_str());
     } else if (name_.find("Vector") != string::npos) {
       TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_DIR3F, adr_,
                       attr_.c_str());
     }
   }
   void send(int program) {
     glUniform3fv(glGetUniformLocation(program, name_.c_str()), 1, adr_);
   }
   bool link(KeyFrame* kf) {
    adr_ = (float*)kf->map_address(type_, name_, 1);
    return adr_ != NULL;
   }
   bool ok() { return adr_ != NULL; }
 private:
     Vec3Uniform(const Vec3Uniform& other) :
             adr_(other.adr_), name_(other.name_),
             type_(other.type_), attr_(other.attr_) {}
   Vec3Uniform& operator=(const Vec3Uniform& other);

   float* adr_;
   string name_;
   string type_;
   string attr_;
};

bool Uniforms::parseLine(const string& line, iUniformPtr* uni) {
  if (line.empty() || line[0] != 'u') return false;

  istringstream is(line);

  string uniform;
  is >> uniform;

  if (uniform.compare("uniform")) return false;

  string type;
  is >> type;

  if (type.compare("int") == 0) {
    iUniformPtr tmp(new IntUniform);
    if (!tmp->parse(line)) return false;
    *uni = tmp;
  } else if (type.compare("float") == 0) {
    iUniformPtr tmp(new FloatUniform);
    if (!tmp->parse(line)) return false;
    *uni = tmp;
#if defined(GL_ARB_gpu_shader_fp64)
  } else if (type.compare("double") == 0) {
    iUniformPtr tmp(new DoubleUniform);
    if (!tmp->parse(line)) return false;
    *uni = tmp;
#endif
  } else if (type.compare("vec3") == 0) {
    iUniformPtr tmp(new Vec3Uniform);
    if (!tmp->parse(line)) return false;
    *uni = tmp;
  } else if (type.compare("bool") == 0) {
    iUniformPtr tmp(new BoolUniform);
    if (!tmp->parse(line)) return false;
    *uni = tmp;
  } else {
    return false;
  }

  return true;
}

bool Uniforms::parseFromGlsl(const string& glsl) {
  istringstream in(glsl);
  string line;
  while (getline(in, line)) {
    iUniformPtr uni;
    if (parseLine(line, &uni)) {
        cout << "glsl UNI: " << uni->toString() << endl;
        uniforms.insert(make_pair(uni->name(), uni));
    }
  }
  return true;
}

void Uniforms::link(KeyFrame* kf) {
  for (unordered_map<string, iUniformPtr>::iterator it =
    uniforms.begin(); it != uniforms.end(); ++it) {
      if (it->second->link(kf)) {
        cout << "link UNI: " << it->second->toString() << endl;
      }
  }
}

void Uniforms::bindToUI(void* bar) {
  for (unordered_map<string, iUniformPtr>::iterator it =
    uniforms.begin(); it != uniforms.end(); ++it) {
      if (it->second->ok())
        it->second->bindToUI(bar);
  }
}

void Uniforms::send(int program) {
  for (unordered_map<string, iUniformPtr>::iterator it =
    uniforms.begin(); it != uniforms.end(); ++it) {
      if (it->second->ok())
        it->second->send(program);
  }
}
