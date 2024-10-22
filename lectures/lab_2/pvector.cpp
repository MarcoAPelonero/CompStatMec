#include "./pvector.hpp"
#include <iostream>

int main(void)
{
  // Il fatto che le parentesi non ci devono stare e perche cpp permette di fare questa cosa, 
  // ma sintatticamente e scorretto

  //pvector vec1{1.0,1.0,2.3}, vec2{2.1,4.0,5.2}, vec3;
  pvector vec1, vec2, vec3;

  for (int i=0; i < NT; i++)
    {  
      vec1.set(i, 0.1*(i+1));
      vec2.set(i, 0.3*(i+1));
    }
  vec1.show("vec1=");
  vec2.show("vec2=");
  vec3 = vec1+vec2;
  vec3.show("vec3=vec1+vec2=");

  (vec3=vec2).set(0,-3);
  vec3.show("vec3=");
  // Oppure posso sovrascrivere 
  vec3 = pvector({0.5, 0.5, 0.5});
  vec3.show("vec3=");
  // O Addirittura 
  vec3 = {0.5, 0.5, 0.5};
  vec3.show("vec3=");

  // Vediamo se funziona
  vec3 = {0.5, 0.5, 0.5, 0.4};
  vec3.show("vec3=");
  // No perche la dimensione dell'array che sta alla base di pvector, è 3, se vado oltre semplicemente non lo salva
  pvector vec4 = {0.5, 0.5, 0.5};
  vec4 += 0.4;
  std::cout << "Exercise" << std::endl;
  vec4.show("vec4=");
  return 0;
}
