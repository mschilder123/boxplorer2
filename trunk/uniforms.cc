#include "uniforms.h"

#include <iostream>
#include <string>
#include <sstream>

#include <AntTweakBar.h>
#include "camera.h"

#define GL_DECLARE_ONLY
#include "shader_procs.h"

using namespace std;

namespace {
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
	 IntUniform(const string& line, KeyFrame* kf) : adr_(NULL) {
		 if (!parseLine(line, &type_, &name_, &attr_)) return;
		 adr_ = (int*)kf->map_address(type_, name_, 1);
	 }
	 virtual ~IntUniform() {}
	 iUniform* Clone() { return new IntUniform(*this); }

 	 const string& name() { return name_; }
	 string toString() {
		ostringstream o;
		o << "int " << name_ << " ";
		if (adr_ != NULL) o << *adr_  << " @" << adr_;
		else o << "(nil)";
		return o.str();
	 }
	 void twVar(void* vbar) {
		 TwAddVarRW((TwBar*)vbar, name_.c_str(), TW_TYPE_INT32, adr_, attr_.c_str());
	 }
	 void send(int program) {
		 glUniform1i(glGetUniformLocation(program, name_.c_str()), *adr_);
	 }
	 bool ok() { return adr_ != NULL; }
 private:
  IntUniform(const IntUniform& other) : adr_(other.adr_), name_(other.name_), type_(other.type_), attr_(other.attr_) {}
  IntUniform& operator=(const IntUniform& other);

  int* adr_;
  string name_;
  string type_;
  string attr_;
};

class FloatUniform : public iUniform {
 public:
	 FloatUniform(const string& line, KeyFrame* kf) : adr_(NULL) {
		if (!parseLine(line, &type_, &name_, &attr_)) return;
		adr_ = (float*)kf->map_address(type_, name_, 1);
	 }
	 virtual ~FloatUniform() {}
	 iUniform* Clone() { return new FloatUniform(*this); }

	 const string& name() { return name_; }
	 string toString() {
		ostringstream o;
		o << "float " << name_ << " ";
		if (adr_ != NULL) o << *adr_ << " @" << adr_;
		else o << "(nil)";
		return o.str();
	 }
	 void twVar(void* bar) {
		 TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_FLOAT, adr_, attr_.c_str());
	 }
	 void send(int program) {
		 glUniform1f(glGetUniformLocation(program, name_.c_str()), *adr_);
	 }
	 bool ok() { return adr_ != NULL; }
private:
    FloatUniform(const FloatUniform& other) : adr_(other.adr_), name_(other.name_), type_(other.type_), attr_(other.attr_) {}
	FloatUniform& operator=(const FloatUniform& other);

	float* adr_;
	string name_;
	string type_;
	string attr_;
};

#if defined(GL_ARB_gpu_shader_fp64)
class DoubleUniform : public iUniform {
 public:
	 DoubleUniform(const string& line, KeyFrame* kf) : adr_(NULL) {
		if (!parseLine(line, &type_, &name_, &attr_)) return;
		adr_ = (double*)kf->map_address(type_, name_, 1);
	 }
	 virtual ~DoubleUniform() {}
	 iUniform* Clone() { return new DoubleUniform(*this); }

	 const string& name() { return name_; }
	 string toString() {
		ostringstream o;
		o << "double " << name_ << " ";
		if (adr_ != NULL) o << *adr_ << " @" << adr_;
		else o << "(nil)";
		return o.str();
	 }
	 void twVar(void* bar) {
		 TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_DOUBLE, adr_, attr_.c_str());
	 }
	 void send(int program) {
		 glUniform1d(glGetUniformLocation(program, name_.c_str()), *adr_);
	 }
	 bool ok() { return adr_ != NULL; }
private:
    DoubleUniform(const DoubleUniform& other) : adr_(other.adr_), name_(other.name_), type_(other.type_), attr_(other.attr_) {}
	DoubleUniform& operator=(const DoubleUniform& other);

	double* adr_;
	string name_;
	string type_;
	string attr_;
};
#endif

class Vec3Uniform : public iUniform {
 public:
	 Vec3Uniform(const string& line, KeyFrame* kf) : adr_(NULL) {
		if (!parseLine(line, &type_, &name_, &attr_)) return;
		adr_ = (float*)kf->map_address(type_, name_, 3);
	 }
	 virtual ~Vec3Uniform() {}
	 iUniform* Clone() { return new Vec3Uniform(*this); }

	 const string& name() { return name_; }
	 string toString() {
		ostringstream o;
		o << "vec3 " << name_ << " ";
		if (adr_ != NULL) o << adr_[0] << " " << adr_[1] << " " << adr_[2] << " @" << adr_;
		else o << "(nil)";
		return o.str();
	 }
	 void twVar(void* bar) {
		 if (name_.find("Color") != string::npos) {
		   TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_COLOR3F, adr_, attr_.c_str());
		 } else if (name_.find("Vector") != string::npos) {
		   TwAddVarRW((TwBar*)bar, name_.c_str(), TW_TYPE_DIR3F, adr_, attr_.c_str());
		 }
	 }
	 void send(int program) {
		 glUniform3fv(glGetUniformLocation(program, name_.c_str()), 1, adr_);
	 }
	 bool ok() { return adr_ != NULL; }
 private:
     Vec3Uniform(const Vec3Uniform& other) : adr_(other.adr_), name_(other.name_), type_(other.type_), attr_(other.attr_) {}
	 Vec3Uniform& operator=(const Vec3Uniform& other);

	 float* adr_;
	 string name_;
	 string type_;
	 string attr_;
};

iUniformPtr link_uniform(const string& line, KeyFrame* kf) {
	if (line.empty() || line[0] != 'u') return iUniformPtr(NULL);

	istringstream is(line);

	string uniform;
	is >> uniform;

	if (uniform.compare("uniform")) return iUniformPtr(NULL);

	string type;
	is >> type;

	if (type.compare("int") == 0) {
		return iUniformPtr(new IntUniform(line, kf));
	} else if (type.compare("float") == 0) {
		return iUniformPtr(new FloatUniform(line, kf));
#if defined(GL_ARB_gpu_shader_fp64)
	} else if (type.compare("double") == 0) {
		return iUniformPtr(new DoubleUniform(line, kf));
#endif
	} else if (type.compare("vec3") == 0) {
		return iUniformPtr(new Vec3Uniform(line, kf));
	}

	return iUniformPtr(NULL);
}
