// from http://shankel.best.vwh.net/interpolate.c

#include <assert.h>
#include <math.h>
#include "interpolate.h"

// Compute the dot product of two vectors.
double dot(const double x[3], const double y[3]) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// Normalize a vector. If it was zero, return false.
bool normalize(double x[3]) {
  double len = dot(x, x); if (len == 0) return false;
  len = 1/sqrt(len); x[0] *= len; x[1] *= len; x[2] *= len;
  return true;
}

void mslerp(const double *m1,const double *m2,double *mr,double t) {
	double q1[4],q2[4],qr[4];
	mat2quat(m1,q1);
	mat2quat(m2,q2);
	qslerp(q1,q2,qr,t);
	quat2mat(qr,mr);
}

void quat2mat(const double *q,double *m) {
	double x2,y2,z2,xy,xz,xw,yz,yw,zw;
	x2 = 2*q[0]*q[0];
	y2 = 2*q[1]*q[1];
	z2 = 2*q[2]*q[2];
	xy = 2*q[0]*q[1];
	xz = 2*q[0]*q[2];
	xw = 2*q[0]*q[3];
	yz = 2*q[1]*q[2];
	yw = 2*q[1]*q[3];
	zw = 2*q[2]*q[3];

	m[0] = 1-y2-z2;
	m[4] = xy+zw;
	m[8] = xz-yw;
//	m[12] = 0;
	
	m[1] = xy-zw;
	m[5] = 1-x2-z2;
	m[9] = yz+xw;
//	m[13] = 0;

	m[2] = xz+yw;
	m[6] = yz-xw;
	m[10] = 1-x2-y2;
//	m[14] = 0;

	m[3] = 0;
	m[7] = 0;
	m[11] = 0;
//	m[15] = 1;
}

void mat2quat(const double *m1, double *q) {
	int i,j,k;
	double Tr;
	double S;
	double m[16];

	// From column major to row major. Only relevant entries are copied.
	m[0] = m1[0];
	m[1] = m1[4];
	m[2] = m1[8];
//	m[3] = 0;
	m[4] = m1[1];
	m[5] = m1[5];
	m[6] = m1[9];
//	m[7] = 0;
	m[8] = m1[2];
	m[9] = m1[6];
	m[10] = m1[10];
//	m[11] = 0;
//	m[12] = m1[12];
//	m[13] = m1[13];
//	m[14] = m1[14];
//	m[15] = m1[15];

	Tr = m[0] + m[5] + m[10] + 1;
	
	if (Tr >= 1) {
		S = 2.0 * sqrt(Tr);
		q[3] = S / 4.0f;
		q[0] = (m[6]-m[9]) / S;
		q[1] = (m[8]-m[2]) / S;
		q[2] = (m[1]-m[4]) / S;
	} else {
		i=0;
		if (m[5] > m[0]) {
			i=1;
		}
		if (m[10] > m[5*i]) {
			i=2;
		}
		j = (i+1)%3;
		k = (j+1)%3;
		
		S = 2.0 * sqrt(1 + m[5*i] - m[5*j] - m[5*k]);
		q[i] = S / 4;
		q[j] = (m[i*4+j] + m[i+4*j]) / S;
		q[k] = (m[i*4+k] + m[i+4*k]) / S;
		q[3] = (m[k+4*j] - m[j+4*k]) / S;
	}
}

void qmul(double* q1, const double* q2) {
	double t[4];

#define x 0
#define y 1
#define z 2
#define w 3

	t[x] = q1[w]*q2[x] + q1[x]*q2[w] + q1[y]*q2[z] - q1[z]*q2[y];
	t[y] = q1[w]*q2[y] + q1[y]*q2[w] + q1[z]*q2[x] - q1[x]*q2[z];
	t[z] = q1[w]*q2[z] + q1[z]*q2[w] + q1[x]*q2[y] - q1[y]*q2[x];
	t[w] = q1[w]*q2[w] - q1[x]*q2[x] - q1[y]*q2[y] - q1[z]*q2[z];

#undef w
#undef z
#undef y
#undef x

	q1[0] = t[0];
	q1[1] = t[1];
	q1[2] = t[2];
	q1[3] = t[3];
}

void qslerp(const double *q1,const double *q2,double *qr,double t) {
	double q3[4];
	int i;
	double sina,sinat,sinaomt,theta;
	double dot = q1[0]*q2[0]+q1[1]*q2[1]+q1[2]*q2[2]+q1[3]*q2[3];
	
	if (dot < 0) {
		dot = -dot;
		for (i=0;i<4;i++) {
			q3[i] = -q2[i];
		}
	}	else {
		for (i=0;i<4;i++) {
			q3[i] = q2[i];
		}
	}
	if (dot > 0.999f) {
        // Very close, just linear interpolate.
		for (i=0;i<4;i++) {
			qr[i] = q1[i] + (q3[i]-q1[i])*t;
		}
	}	else {
		theta = acos(dot);
		sina = sin(theta);
		sinat = sin(theta*t);
		sinaomt = sin(theta*(1-t));
		for (i=0;i<4;i++) {
			qr[i] = (q1[i]*sinaomt+q3[i]*sinat)/sina;
		}
	}
}

void qnormalize(double* q) {
  double mag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
//  if (fabs(mag) > .000001 && fabs(mag - 1.0) > .000001)
  {
    double invMag = 1.0 / mag;
    q[0] *= invMag;
    q[1] *= invMag;
    q[2] *= invMag;
    q[3] *= invMag;
  }
}

void qinvert(double* qout, const double* q) {
  qout[0] = -q[0];
  qout[1] = -q[1];
  qout[2] = -q[2];
  qout[3] = q[3];
  qnormalize(qout);
}

// from game programming gems 2, page 224
void quat2x(const double* q, double* x) {
  double d = sqrt(2*(1-q[3]));
  // TODO: deal with d ~ 0
  x[0] = q[0] / d;
  x[1] = q[1] / d;
  x[2] = q[2] / d;
  x[3] = (1-q[3]) / d;
}

void x2quat(const double* x, double *q) {
  double d = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
  // TODO: deal with d ~ 0
  q[3] = (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - x[3]*x[3]) / d;
  q[0] = (2*x[0]*x[3]) / d;
  q[1] = (2*x[1]*x[3]) / d;
  q[2] = (2*x[2]*x[3]) / d;
}
