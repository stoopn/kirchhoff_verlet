#ifndef QUATERNION_H
#define QUATERNION_H 1

#include "matrix.h"


class Quaternion
{
public:
  double q[4];
//  int t;
  Quaternion();
  Quaternion(double q0, double q1, double q2, double q3);
  Quaternion conj(void);
  Quaternion operator * (Quaternion);
  Quaternion operator + (Quaternion);
  Quaternion operator - (Quaternion);
  Quaternion operator * (double);
  double dot(Quaternion);


  void set(double q0, double q1, double q2, double q3);
  void set(double q3, Tvect v);
  void set(double phi, double theta, double psi);
  void set(float ql[4]);
  void reset();

  Tvect getD3();
  Tvect getdD3dq0();
  Tvect getdD3dq1();
  Tvect getdD3dq2();
  Tvect getdD3dq3();

  void buildRmat(Tmat& R, int transpose);

  Quaternion rotate(Tmat rotM);
  double norm();
  void normalize();
  void display();
  
  Tvect transform(Tvect v, int t);
};

//Quaternion acelQ(Quaternion Q, Quaternion QV, Tvect wa);
Tvect acelT(Quaternion Q, Quaternion t);

#endif
