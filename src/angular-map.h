#pragma once

#include "flattened_mesh.h"

namespace flattener
{

class ANGULAR_MAP_CREATOR
{
public:
  virtual bool create(masc::model * m, flattened_mesh& fm) = 0;

protected:
  virtual bool cut(masc::model * m);
  bool cut(masc::model * m, uint seed);
};

class CGAL_AMAP_Creator : public ANGULAR_MAP_CREATOR
{
public:

  enum Method {Tutte, Conformal_Map, Floater, Authalic, LeastSquares};
  enum Boundary_Type {Circular_Arc_Length, Circular_Uniform, Square_Arc_Length, Square_Uniform};

  CGAL_AMAP_Creator(Method method, Boundary_Type bd):param_method_(method), bd_type_(bd){}
  virtual ~CGAL_AMAP_Creator(){}

  virtual bool create(masc::model * m, flattened_mesh& fm);

private:

  Method param_method_;
  Boundary_Type bd_type_;
};


class LP_AMAP_Creator : public ANGULAR_MAP_CREATOR
{
public:

  LP_AMAP_Creator(){}
  virtual ~LP_AMAP_Creator(){}

  virtual bool create(masc::model * m, flattened_mesh& fm);

protected:

};


} //end namespace flattener
