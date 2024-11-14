#include<iostream>
#include "./randnumgen.hpp"
int main(void)
{
  rng.rseed(); // random seed, l'equivalente di srand48(time(NULL))
  //rng.seed(1); // fixed seed
  std::cout << "random numbers in [0,1):\n";
  for (int i=0; i < 10; i++)
    {
      std::cout << rng.ranf();
      if (i < 9) 
        std::cout << ", ";
      else 
        std::cout << "\n";
    }
}
