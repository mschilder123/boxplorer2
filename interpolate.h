#ifndef _F_INTERPOLATE_H_
#define _F_INTERPOLATE_H_

#ifdef PI
  #undef PI
#endif
#define PI          3.14159265358979324

// Compute the dot product of two vectors.
double dot(const double x[3], const double y[3]);

// Normalize a vector. If it was zero, return false.
bool normalize(double x[3]);

// SLERP 2 OpenGL matrices
void mslerp(const double *m1,const double *m2,double *mr,double t);

// SLERP 2 quaternions
void qslerp(const double *q1,const double *q2,double *qr,double t);

// Matrix-Quaternion conversions
void mat2quat(const double *m,double *q);
// Does not touch m[12..14]
void quat2mat(const double *q,double *m);

// Quat to splinable r4 mapping conversions.
// Spline over x to get c2 continuous quat spline.
void quat2x(const double* q, double* x);
void x2quat(const double* x, double *q);

void qnormalize(double* q);
void qinvert(double* qout, const double* q);

// q1 *= q2
void qmul(double* q1, const double* q2);

#endif
