#ifndef MATRIX_H__
#define MATRIX_H__

#include <cassert>
#include <memory.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <memory.h>

using std::ostream;
using std::endl;
using std::cout;
using std::vector;

/**
 * matrix library for CS171, lab 1
 */

/*** forward declarations ***/
template <int M, int N>
class Matrix;

/*** matrix declarations ***/
typedef Matrix<4, 4> Matrix44;
typedef Matrix<3, 3> Matrix33;

struct MatrixId { };

template <int M, int N>
class Matrix
{
protected:
  double mat[M*N];

  // unsafe element access
  double elem_c(int i, int j) const { return mat[i * N + j]; }
  double& elem(int i, int j) { return mat[i * N + j]; }

  /*** vectors only ***/
  double elem_c(int i) const { return mat[i]; }

public:
  Matrix<M, N>() { } 
  Matrix<M, N>(double s) { fill(s); }
  Matrix<M, N>(MatrixId ic) { identity(); }
  Matrix<M, N>(double* s) { memcpy(mat, s, sizeof(double)*M*N); }

  // assignment operator: default is acceptable
  // copy constructor: default is acceptable
  // destructor: defualt is acceptable

  // construction shortcuts
  // fills the matrix with zeros
  void fill(double s)
  {
    for (int i = 0; i < M*N; ++i)
      mat[i] = s;
  }
  void zeros() { fill(0); }
  
  // generalization of the identity to non-square matrices
  // one's along the main diagonal, zeroes elsewhere 
  void identity() 
  {
    for (int i = 0; i < M; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        if (i == j) 
          elem(i, j) = 1;
        else
          elem(i, j) = 0;
      }
    }
  }

  /***************************************
   * indexing
   **************************************/
  
  // operator() returns references, and so are valid lhs arguments
  // see operator[] for const-safe alternative
  double& operator()(int i, int j)
  {
    assert(0 <= i && i < M);
    assert(0 <= j && j < N);
    return elem(i, j);
  }

  double operator()(int i, int j) const
  {
    assert(0 <= i && i < M);
    assert(0 <= j && j < N);
    return elem_c(i, j);
  }

  double& operator()(int i)
  {
    //assert(N == 1);
    assert(0 <= i && i < M);
    return elem(i, 0);    
  }

  double operator()(int i) const
  {
    assert(0 <= i && i < M);
    return elem_c(i, 0);
  }

  /***************************************
   * simple arithmetic functions
   **************************************/
  // addition
  Matrix<M, N>& operator+=(const Matrix<M, N>& m)
  {
    for (int i = 0; i < M*N; ++i)
      mat[i] += m.mat[i];
  
    return *this;
  }

  Matrix<M, N> operator+(const Matrix<M, N>& m) const
  {
    Matrix<M, N> ret = *this;
    return (ret += m);
  }
  
  // subtraction
  Matrix<M, N>& operator-=(const Matrix<M, N>& m)
  {
    for (int i = 0; i < M*N; ++i)
      mat[i] -= m.mat[i];
  
    return *this;
  }

  Matrix<M, N> operator-(const Matrix<M, N>& m) const
  {
    Matrix<M, N> ret = *this;
    return (ret -= m);
  }

  // negation
  Matrix<M, N> operator-() const
  {
    return (*this) * -1.0f;
  }

  // double multiplication
  Matrix<M, N>& operator*=(double s)
  {
    for (int i = 0; i < M*N; ++i)
      mat[i] *= s;
  
    return *this;
  }

  template <int A, int B>
  friend Matrix<A, B>& operator*=(double s, Matrix<A, B>& m);

  Matrix<M, N> operator*(double s) const
  {
    Matrix<M, N> ret = *this;
    return (ret *= s);
  }
  
  template <int A, int B>
  friend Matrix<A, B> operator*(double s, Matrix<A, B> m);

  // double division
  Matrix<M, N>& operator/=(double s)
  {
    for (int i = 0; i < M*N; ++i)
      mat[i] /= s;
  
    return *this;
  }

  Matrix<M, N> operator/(double s) const
  {
    Matrix<M, N> ret = *this;
    return (ret /= s);
  }


  // matrix multiplication
  template<int A, int B, int C>
  friend Matrix<A, C> operator*(const Matrix<A, B>& m1, const Matrix<B, C>& m2);
  
  template<int P>
  Matrix<M, P>& operator*=(const Matrix<N, P>& m)
  {
    return (*this = *this * m);
  }

  // row operations
  void scale(int row, double s)
  {
    assert(0 <= row);
    assert(row < M);
    for (int j = 0; j < N; ++j)
      elem(row, j) *= s;
  }

  void swap(int r1, int r2)
  {
    assert(0 <= r1 && r1 < M);
    assert(0 <= r2 && r2 < M);
    
    if (r1 == r2)
      return;

    double temp;
    for (int j = 0; j < N; ++j)
    {
      temp = elem_c(r1, j);
      elem(r1, j) = elem_c(r2, j);
      elem(r2, j) = temp;
    }
  }
  
  void mrowadd(int r1, int r2, double s)
  {
    assert(0 <= r1 && r1 < M);
    assert(0 <= r2 && r2 < M);
    
    for (int j = 0; j < N; ++j)
    {
      elem(r1, j) += elem_c(r2, j) * s;
    }
  }
  
  /*** SQAURE MATRIX OPERATIONS ***/
  
  // inversion
  // only valid for square matrices
  void invert()
  {
    assert(M == N);
    Matrix<M, N> res;
    res.identity();

    for (int r = 0; r < M; ++r)
    {
      // find the row with the largest pivot and swap it to the top
      // if the current pivot is 0
      int pivot_row = r;

      for (int s = r + 1; s < M; ++s)
      {
        if (abs(elem_c(pivot_row, r)) < abs(elem_c(s, r)))
          pivot_row = s;
      }

      this->swap(pivot_row, r);
      res.swap(pivot_row, r);
 
      double pivot = elem_c(r, r);
      
      // scale the row, and perform the appropriate mult-row-adds
      scale(r, 1 / pivot);
      res.scale(r, 1 / pivot);

      for (int s = 0; s < M; ++s)
      {
        if (s != r)
        {
          double f = -elem_c(s, r);
          mrowadd(s, r, f);
          res.mrowadd(s, r, f);
        }
      }
    }

    // copy the resulting matrix back into the original
    *this = res;
  }

  Matrix<M, N> inverse() const
  {
    Matrix<M, N> inv = *this;
    inv.invert();
    return inv;
  }

  /***************************************
   * resizing operations 
   **************************************/
  /*** return the same matrix extended by one row and one column,
       with 0's on the side and 1 in the corner */
  Matrix<M+1, N+1> matrix_extend() const
  {
    Matrix<M+1, N+1> res;

    for (int i = 0; i < M; ++ i)
    {
      for (int  j = 0; j < N; ++j)
      {
        res(i, j) = elem_c(i, j);
      }
      res(i, N) = 0.0;
    }

    for (int j = 0; j < N; ++j)
      res(M, j) = 0.0;

    res(M, N) = 1.0;
    return res;
  }

  /*** extend a column vector by one specified element ***/
  Matrix<M+1, 1> vector_extend(double f = 0) const
  {
    assert(N == 1);
    Matrix<M+1, 1> res;
    for (int i = 0; i < M; ++i)
      res(i, 0)  = elem_c(i, 0);
    res(M, 0) = f;
   
    return res;
  }

  Matrix<M-1, N-1> matrix_chop() const
  {
    Matrix<M-1, N-1> res;
    for (int i = 0; i < M-1; ++i)
    {
      for (int j = 0; j < N-1; ++j)
      {
        res(i, j) = elem_c(i, j);
      }
    }

    return res;
  }        

  Matrix<M-1, 1> vector_chop() const
  {
    Matrix<M-1, 1> res;
    for (int i = 0; i < M-1; ++i)
    {
      res(i, 0) = mat[i];
    }

    return res;
  }

  /**************************************
   * vector operations
   *************************************/
  double dot(const Matrix<M, 1>& v) const
  {
    double result = 0;
    for (int i = 0; i < M; ++i)
      result += mat[i] * v.mat[i];
    
    return result;
  }

  // specifically for 3-vectors
  Matrix<3, 1> cross(const Matrix<3, 1>& v) const
  {
    assert(M == 3);
    assert(N == 1);
    Matrix<3, 1> res;
    res.mat[0] = mat[1] * v.mat[2] - mat[2] * v.mat[1];
    res.mat[1] = mat[2] * v.mat[0] - mat[0] * v.mat[2];
    res.mat[2] = mat[0] * v.mat[1] - mat[1] * v.mat[0];

    return res;
  }

  // return [a], where [a] . b = a x b
  Matrix<3, 3> cross_gradient() const
  {
    double res[9] = { 0,      -mat[2], mat[1],
                      mat[2],  0,     -mat[0],
                      -mat[1], mat[0],  0      };
    
    return Matrix<3, 3>(res);
  }

  // shorthand for outer product with another vector
  Matrix<M, M> outer(const Matrix<M, 1>& v) const
  {
    return *this * v.transpose();
  }

  Matrix<M, 1>& homogenize()
  {
    assert(N == 1);
    double s = mat[M - 1];
    return (*this /= s);
  }

  double mag_sq() const
  {
    return dot(*this);
  }

  double mag() const
  {
    return sqrt(mag_sq());
  }
  
  Matrix<M, 1> normal() const
  {
    return *this / mag();
  }

  Matrix<M, 1>& normalize()
  {
    return (*this = this->normal());
  }

  // operator[] returns values, and are const-safe
  double operator[](int i) const
  {
    return this->operator()(i);
  }


  /**************************************
   * miscellaneous operations 
   *************************************/
  Matrix<N, M> transpose() const
  {
    Matrix<N, M> ret;
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        ret(j, i) = (*this)(i, j);

    return ret;
  }

  static Matrix<M, N> id;
  static Matrix<M, N> zero;
};

/***************************************
 * friend operators 
 **************************************/
template <int A, int B>
Matrix<A, B>& operator*=(double s, Matrix<A, B>& m)
{
  return m *= s;
}

template <int A, int B>
Matrix<A, B> operator*(double s, Matrix<A, B> m)
{
  return m * s;
}

template<int M, int N, int P>
Matrix<M, P> operator*(const Matrix<M, N>& m1, const Matrix<N, P>& m2)
{
  Matrix<M, P> res;
  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < P; ++j)
    {
      res(i, j) = 0;
      for (int k = 0; k < N; ++k)
      {
        res(i, j) += m1.elem_c(i, k) * m2.elem_c(k, j);
      }
    }
  }
  
  return res;
}

/***************************************
 * output functions
 **************************************/

template <int A, int B>
ostream& operator <<(ostream& o, const Matrix<A, B>& m)
{
  for (int i = 0; i < A; ++i)
  {
    for (int j = 0; j < B; ++j)
    {
      o.width(10);
      o.fill(' ');
      o.precision(5);
      o << m(i, j);
    }
    o << endl;
  }
  return o;
}

/** vectors are output specially, by row **/
template <int M>
ostream& operator <<(ostream& o, const Matrix<M, 1>& m)
{

  o << '{' << m(0);
  for (int i = 1; i < M; ++i)
    o << ", " << m(i);
  o << '}';
  return o;
}

/***************************************
 * vector functions
 **************************************/
typedef Matrix<2, 1> Vector2;
typedef Matrix<3, 1> Vector3;
typedef Matrix<4, 1> Vector4;
typedef Vector3 Tvect;
typedef Matrix33 Tmat;

// vector construction shortcuts
Vector2 Vec(double x, double y);
Vector3 Vec(double x, double y, double z);
Vector4 Vec(double x, double y, double z, double w);

#endif

