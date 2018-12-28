#include "quaternion.h"

Quaternion::Quaternion()
{
  int k;
  for( k=0; k< 4; k++) q[k]=0.0;
}

Quaternion::Quaternion(double q0, double q1, double q2, double q3)
{
  q[0] = q0;
  q[1] = q1;
  q[2] = q2;
  q[3] = q3;
}
void Quaternion::reset()
{
  q[0] = 0.;
  q[1] = 0.;
  q[2] = 0.;
  q[3] = 0.;
}

void Quaternion::set(double q0, double q1, double q2, double q3)
{
  q[0] = q0;
  q[1] = q1;
  q[2] = q2;
  q[3] = q3;
}

void Quaternion::set(double phi, double theta, double psi)
{
  double a1, a2, a3;

  a1 = 0.5*theta;
  a2 = 0.5*(phi-psi);
  a3 = 0.5*(phi+psi);

  q[1] = sin(a1)*cos(a2);
  q[2] = sin(a1)*sin(a2);
  q[3] = cos(a1)*sin(a3);
  q[0] = cos(a1)*cos(a3);
}

void Quaternion::set(float q1[4])
{
	q[0] = (double) q1[0];
	q[1] = (double) q1[1];
	q[2] = (double) q1[2];
	q[3] = (double) q1[3];
}

void Quaternion::set(double q3, Tvect v)
{
  for(int k = 1; k<4; k++) q[k] = v[k];
  q[0] = q3;
}

double Quaternion::norm()
{	
  return (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void Quaternion::normalize()
{
    double len = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if(len>0.0) 
	for(int k=0;k<4;k++) q[k] = q[k]/len;
}

Quaternion Quaternion::operator * (Quaternion q1)
{
  Quaternion rq;

  rq.q[0] = q[0]*q1.q[0] - q[1]*q1.q[1] - q[2]*q1.q[2] - q[3]*q1.q[3];
  rq.q[1] = q[1]*q1.q[0] + q[0]*q1.q[1] - q[3]*q1.q[2] + q[2]*q1.q[3];
  rq.q[2] = q[2]*q1.q[0] + q[0]*q1.q[2] - q[1]*q1.q[3] + q[3]*q1.q[1];
  rq.q[3] = q[3]*q1.q[0] + q[0]*q1.q[3] + q[1]*q1.q[2] - q[2]*q1.q[1];

  return rq;
}

double Quaternion::dot(Quaternion q1)
{
  return q[0]*q1.q[0]+ q[1]*q1.q[1]+ q[2]*q1.q[2]+ q[3]*q1.q[3];
}

Tvect Quaternion::getD3()
{
  Tvect d3;
  d3(0) = 2.*(q[0]*q[2] + q[1]*q[3]);
  d3(1) = 2.*(q[1]*q[2] - q[0]*q[3]);
  d3(2) = -q[0]*q[0] - q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
  return d3;
}

Tvect Quaternion::getdD3dq0()
{
  Tvect dd;
  dd(0) = 2.*q[2];
  dd(1) = -2.*q[3];
  dd(2) = -2.*q[0];
  return dd;
}

Tvect Quaternion::getdD3dq1()
{
  Tvect dd;
  dd(0) = 2.*q[3];
  dd(1) = 2.*q[2];
  dd(2) = -2.*q[1];
  return dd;
}

Tvect Quaternion::getdD3dq2()
{
  Tvect dd;
  dd(0) = 2.*q[0];
  dd(1) = 2.*q[1];
  dd(2) = 2.*q[2];
  return dd;
}

Tvect Quaternion::getdD3dq3()
{
  Tvect dd;
  dd(0) = 2.*q[1];
  dd(1) = -2.*q[0];
  dd(2) = 2.*q[3];
  return dd;
}

Quaternion Quaternion::operator * (double scal)
{
  Quaternion rq;
	
  for(int k = 0; k<4; k++) rq.q[k] = q[k]*scal;
  

  return rq;
}


Quaternion Quaternion::operator + (Quaternion q1)
{
  Quaternion rq;
  for(int k = 0; k<4; k++) rq.q[k] = q[k] + q1.q[k];
 
  return rq;
}

Quaternion Quaternion::operator - (Quaternion q1)
{
  Quaternion rq;
  for(int k = 0; k<4; k++) rq.q[k] = q[k] - q1.q[k];
 
  return rq;
}

Quaternion Quaternion::conj(void)
{
  Quaternion rq;
  for(int k = 1; k<4; k++) rq.q[k] = -q[k];
  rq.q[0] = q[0];
  
  return rq;
}

void Quaternion::display()
{
  std::cout << q[0] << ", " << q[1] << ", " <<q[2] << ", " <<q[3] << "\n";
}

Tvect acelT(Quaternion Q, Quaternion t)
{  
  Quaternion tt = Q.conj()*t*0.5;
  Tvect v;
  v(0) = tt.q[1];
  v(1) = tt.q[2];
  v(2) = tt.q[3];
  return v;
}
