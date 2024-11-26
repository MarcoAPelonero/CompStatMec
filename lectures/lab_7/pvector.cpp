#include "./pvector.hpp"
int main(void)
{
  pvector<double,3> vec1={1,1.0,2.3}, vec2{2.1,4.0,5.2}, vec3;
  double val=2.2;
  pvector<int,3> iv={1,2,3};


  std::cout << "norm(iv)=" << iv.norm() << "\n";
  // seeding the Mersenne-Twister
  rng.rseed();
 
  vec1.show("vec1=");
  vec2.show("vec2=");
  vec3 = vec1+vec2;
  
  vec3.show("vec3=vec1+vec2=");

  // RETURN VALUE BY REFERENCE IN AN ASSIGNMENT OPERATOR
  // = operator returns object vec3 by reference,
  // hence method set will set first element of vec3 to -3
  (vec3=vec2).set(0,-3);
  vec3.show("vec3=");

  val = vec1(2);
  std::cout << "val=" << val << std::endl;

  char ival=3;
  
  vec3 = {ival,val,1};
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>> vec3=" << vec3 << std::endl;
  // CONST VECTOR WITH operator ()
  // if const operator() is not defined we will obtain
  // an error using operator () with const pvector 
  const pvector<double,3> vv={1.0,2.0,6.3};
  std::cout << "norm(vv)=" << vv.norm() << std::endl;
  vec3(0) = vv(1);
  vec3.show("vec3=");
  vec3 = 1.5*vec2;
  vec3.show(">>>>vec3=");
  
  // element-wise operations
  //
  std::cout << "vec3.mulcw(vec3)" << vec3.mulcw(vec3) << std::endl;
  std::cout << "vec3.divcw(vec3)" << vec3.divcw(vec3) << std::endl;

  // NARROWING
  int a=3;
  // if remove conversion to double of integer a, we get
  // an error here since int can not be "narrowed" to double.
  // This is a peculiar behavior of list initialization
  pvector<double,3> vec4={double(a)+0.3,1.5,22.5};

  vec4.random_orient();
  // overloading of << operator and use of rint
  std::cout << "vec4=" << vec4 <<  " norm(vec4)=" << vec4.norm() << std::endl;
  //
  std::cout << "rint(vec4)=" << rint(vec4) << std::endl;
  return 0;
}
