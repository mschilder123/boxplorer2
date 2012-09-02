#ifdef __cplusplus
extern "C" {
#endif

//SLERP 2 OpenGL matrices
void mslerp(const double *m1,const double *m2,double *mr,double t);

//SLERP 2 quaternions
void qslerp(const double *q1,const double *q2,double *qr,double t);

//Matrix-Quaternion conversions
void mat2quat(const double *m,double *q);
void quat2mat(const double *q,double *m);

//Quat to splinable r4 mapping conversions.
//Spline over x to get c2 continuous quat spline.
void quat2x(const double* q, double* x);
void x2quat(const double* x, double *q);

#ifdef __cplusplus
}
#endif