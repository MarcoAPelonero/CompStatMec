#include<iostream>
#include<chrono>
#include "./randnumgen.hpp"

randnumgen<double,std::mt19937> mt;
randnumgen<double,std::mt19937_64> mt64;
randnumgen<double,std::ranlux24> rl24;
randnumgen<double,std::ranlux48> rl48;

int main(void)
{
  const int NT=5000000;
  using cl=std::chrono::high_resolution_clock;


  // Mersenne-Twister
  //rng.rseed(); // random seed, l'equivalente di srand48(time(NULL))
  mt.seed(1); // fixed seed
  auto tini=cl::now();
  std::cout << "Testing Mersenne-Twister...";
  long int a=0;
  for (int i=0; i < NT; i++)
    {
      a+=10*mt.ranf();
    }
  std::cout << "(a=" << a << ")" << std::endl;
  auto delt=std::chrono::duration_cast<std::chrono::microseconds>(cl::now() - tini).count();
  std::cout << "time elapsed for Mersenne-Twister: " << delt/1.0E6 << " seconds\n\n";
  auto delt0 = delt;

  // Mersenne-Twister64
  mt64.seed(1); // fixed seed
  tini=cl::now();
  std::cout << "Testing Mersenne-Twister64...";
  a=0;
  for (int i=0; i < NT; i++)
    {
      a+=10*mt64.ranf();
    }
  std::cout << "(a=" << a << ")" << std::endl;
  delt=std::chrono::duration_cast<std::chrono::microseconds>(cl::now() - tini).count();
  std::cout << "time elapsed for Mersenne-Twister64: " << delt/1.0E6 << " seconds\n";
  std::cout << "which is " << (double)delt/delt0 << " times slower than Mersenne-Twister\n\n";

  // Ranlux24
  rl24.seed(1); // fixed seed
  tini=cl::now();
  std::cout << "Testing Ranlux24...";
  for (int i=0; i < NT; i++)
    {
      a+=10*rl24.ranf();
    }
  std::cout << "(a=" << a << ")" << std::endl;
  delt=std::chrono::duration_cast<std::chrono::microseconds>(cl::now() - tini).count();
  std::cout << "time elapsed for Ranlux24: " << delt/1.0E6 << " seconds\n";
  std::cout << "which is " << (double)delt/delt0 << " times slower than Mersenne-Twister\n\n";

  // Ranlux48
  rl48.seed(1); // fixed seed
  tini=cl::now();
  std::cout << "Testing Ranlux48...";
  a=0;
  for (int i=0; i < NT; i++)
    {
      a+=10*rl48.ranf();
    }
  std::cout << "(a=" << a << ")" << std::endl;
  delt=std::chrono::duration_cast<std::chrono::microseconds>(cl::now() - tini).count();
  std::cout << "time elapsed for Ranlux48: " << delt/1.0E6 << " seconds\n";
  std::cout << "which is " << (double)delt/delt0 << " times slower than Mersenne-Twister\n";
}
