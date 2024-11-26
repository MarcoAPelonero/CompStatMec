#include"./pmatrix.hpp"
int main(void)
{
  pmatrixq<double,3> M{1,2,3,4,5,6,7,8.3,9}, B, C;
  pvector<double, 3> v={1,2,2};
  B = M;

  C = {2,3,2,3,4,5,8,0,1};
  M = B - C;
  M += B;
  M /= 1.5;
  std::cout << "M=" << M << std::endl;
  M = M.transpose();
  std::cout << "M=" << M << std::endl;
  if (!(M==C))
    std::cout << "diverse\n";
  std::cout << "v=" << v << std::endl;
  std::cout << "M=" << M << std::endl;
  std::cout << (1.2*v)*((M*1.2)*(1.3*M)) << std::endl;
  std::cout << "B=" << B << std::endl;
  std::cout << "C=" << C << std::endl;
  std::cout << "M = B + C = " << M << std::endl;
  return 0;
}
