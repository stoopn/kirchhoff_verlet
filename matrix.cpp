#include "matrix.h"


Vector2 Vec(double x, double y)
{
  Vector2 v;
  v(0) = x; v(1) = y;
  return v;
}


Vector3 Vec(double x, double y, double z)
{
  Vector3 v;
  v(0) = x; v(1) = y; v(2) = z;
  return v;
}

Vector4 Vec(double x, double y, double z, double w)
{
  Vector4 v;
  v(0) = x; v(1) = y; v(2) = z; v(3) = w;
  return v;
}

