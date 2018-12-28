#ifndef UTILS_H
#define UTILS_H 1

#include "matrix.h"
#include "quaternion.h"
#include <deque>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>


#define PI 3.1415927
#define MAX_ANGLE 2.1  // 120 degrees (equal length triangle)
#define MAX_FORCE 200.
#define STICKY_MIN_VEL 0.01
#define RETURN_RGB(r, g, b) {RGB.R = r; RGB.G = g; RGB.B = b;} 
#define KB 1.0




#define CONT_SPHERE 1
namespace blitz {
  inline double dot(Tvect a, Tvect b)
  {
    return a.dot(b);	
  }
  
  inline Tvect cross(Tvect a, Tvect b)
  {
    return a.cross(b);
  }
  inline double pow2(double a)
  {
    return a*a;
  }
  inline double pow3(double a)
  {
    return a*a*a;
  }
  inline double pow4(double a)
  {
    return pow2(a)*pow2(a);
  }
  inline Quaternion mprod(Matrix44 mat, Quaternion q)
  {
    Vector4 v;
    v(0) = q.q[0];
    v(1) = q.q[1];
    v(2) = q.q[2];
    v(3) = q.q[3];
    v = mat*v;
    return Quaternion(v[0],v[1],v[2],v[3]);
  }
  
}

inline Quaternion B1mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(-q.q[3], -q.q[2], q.q[1], q.q[0]);
  else
    return Quaternion(q.q[3], q.q[2], -q.q[1], -q.q[0]);
}
inline Quaternion B2mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(q.q[2], -q.q[3], -q.q[0], q.q[1]);
  else
    return Quaternion(-q.q[2], q.q[3], q.q[0], -q.q[1]);
}
inline Quaternion B3mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(-q.q[1], q.q[0], -q.q[3], q.q[2]);
  else
    return Quaternion(q.q[1], -q.q[0], q.q[3], -q.q[2]);
}

inline Quaternion B01mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(-q.q[3], q.q[2], -q.q[1], q.q[0]);
  else
    return Quaternion(q.q[3], -q.q[2], q.q[1], -q.q[0]);
}
inline Quaternion B02mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(-q.q[2], -q.q[3], q.q[0], q.q[1]);
  else
    return Quaternion(q.q[2], q.q[3], -q.q[0], -q.q[1]);
}
inline Quaternion B03mult(Quaternion q, bool transp=false)
{
  if (transp)
    return Quaternion(q.q[1], -q.q[0], -q.q[3], q.q[2]);
  else
    return Quaternion(-q.q[1], q.q[0], q.q[3], -q.q[2]);
}


/*
inline Tmat cross_gradient(Tvect &mat)
{    
	Tmat mt;
	mt(0,0)=0.;
	mt(0,1)=mat[2];
	mt(0,2)=-mat[1];
	mt(1,0)=-mat[2];
	mt(1,1)=0.;
	mt(1,2)=mat[0];
	mt(2,0)=mat[1];
	mt(2,1)=-mat[0];
	mt(2,2)=0.;
	
	return mt;
}

inline Tmat outer_prod(Tvect &a1, Tvect &a2)
{
	Tmat mt;
	mt(0,0)= a1[0]*a2[0];
	mt(0,1)= a1[0]*a2[1];
	mt(0,2)= a1[0]*a2[2];
	mt(1,0)= a1[1]*a2[0];
	mt(1,1)= a1[1]*a2[1];
	mt(1,2)= a1[1]*a2[2];
	mt(2,0)= a1[2]*a2[0];
	mt(2,1)= a1[2]*a2[1];
	mt(2,2)= a1[2]*a2[2];
	
	return mt;	
}

inline Tvect transform_vect(Tvect v0, Tmat rotM)
{//Rotate the vector, not the coordinate system...

  return blitz::product(rotM,v0);
}
*/

typedef struct {
    float R;
    float G;
    float B;
} RGBType;


inline float finvsqrt(float x) { //Calculates 1/sqrt(x) approximatively
    float xhalf = 0.5f*x;
    int i = *(int*)&x; // get bits for floating value
    i = 0x5f3759df - (i>>1); // give initial guess y0
    x = *(float*)&i; // convert bits back to float
    x *= 1.5f - xhalf*x*x; // newton step, repeating this step
    return x;   
}


inline RGBType HSV_to_RGB(float h, float s, float v) {

    // H is given on [0, 6] or UNDEFINED. S and V are given on [0, 1].
    // RGB are each returned on [0, 1].
    float m, n, f;
    int i;
    RGBType RGB;

    if(h == -1) 
	RETURN_RGB(v, v, v);
    i = floor(h);
    f = h - i;
    if(!(i & 1)) 
	f = 1 - f; // if i is even
    m = v * (1 - s);
    n = v * (1 - s * f);
    switch (i) {
	case 6: RETURN_RGB(v, n, m);
	case 0: RETURN_RGB(v, n, m);
	case 1: RETURN_RGB(n, v, m);
	case 2: RETURN_RGB(m, v, n);
	case 3: RETURN_RGB(m, n, v);
	case 4: RETURN_RGB(n, m, v);
	case 5: RETURN_RGB(v, m, n);
    }
    return RGB;
}


inline double rand01() {
    return ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
}
inline double uniformRandom()
{
    return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}
// return a normally distributed random number
inline double normalRandom(double m, double s)
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1, y2;
    
    
	
    do {
        x1 = 2.0 * uniformRandom() - 1.0;
        x2 = 2.0 * uniformRandom() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    
	return( m + y1 * s );
}

template<typename T>
class MovAvg {
 public:
	T avg;
	
	MovAvg(int size, T aver) {
	data.resize(0);
	MSize=size;
	avg=aver;
	std::cout<<"init1 movav avg: "<<avg<<"\n";	}

    MovAvg(int size) {
	data.resize(0);
	MSize=size;
	avg=0;
	std::cout<<"init2 movav avg: "<<avg<<"\n";	}

    MovAvg() {
	data.resize(0);
	MSize=0;
	avg=0;
	std::cout<<"init3 movav avg: "<<avg<<"\n"; }
    
    
   void operator<<(T value) {
	data.push_back(value);
	if (value<0) {
	    std::cout<<"WARNING: new val: "<<value<<", ";
		std::cout<<"old average: "<<avg<<"\n";
	}
	if (data.size()<=MSize) {
	    avg=(avg*(MSize-1)+value)/MSize;
		std::cout<<"building average: "<<avg<<", size "<<data.size()<<"\n";
	} else if (data.size()>MSize) {
	   avg+=(data.back()-data.front())/MSize;
	   data.pop_front();

	} 
	}
	
    void resize(int newsize) {
	MSize=newsize;
	std::cout<<"resize1 movav avg: "<<avg<<" size: "<<newsize<<"\n"; }
	
    void resize(int newsize, T aver) {
	MSize=newsize;
	avg=aver;
	std::cout<<"resize2 movav avg: "<<avg<<"\n"; }
	
 private:
    int MSize;
    std::deque<T> data;
};





#endif
