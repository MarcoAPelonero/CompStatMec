#ifndef _PVECTOR_
#define _PVECTOR_
// EXERCISE create an efficient vector class in C++ with methods and overloaded operators which you can find in the TODO LIST.
// TODO LIST:
// * = higher priority, X = done, E= exerise
// [X] standard constructor
// [*] constructor (with overloading->initializer_list)
// [X] destructor (empty)
// [*] show()
// [X] overloading of = assignment 
// [X] overloading of + operator for addition of vectors
// [*] overloading of - operator for substraction of vectors
// [*] sum()
// [*] get( ) (get i-th element)
// [*] set( , ) (set i-th element)
// [E] overloading of += and -= operators with vectors
// -------------------- FIN QUI: LAB #1 E #2
//
// [*] overloading of (...) operator (to use vector as C arrays, e.g. v(0)=1
// [*] comma initialization: overloading of << and , operators
// [*] "scalar time vector" through friend function
// [*] overloading of * operator for "vector times scalar" product
// [*] norm() norm of vector, implemented by using overloaded operator "*"
// [*] overloading of operator == (check whether two vectors are equal)
// [*] overloading of *= and /= (with scalars)
// [*] overloading of * operator for scalar product
// [*] overloading of ^ operator for cross product 
// [*] norm(), calculate the modulus of a vector
// --------------------- FIN QUI: LAB #3
//
// [ ] rint 
// [ ] overloading of << for output (friend function)
// [ ] mulcw, multiplication element wise
// [ ] divcw, division element wise
// [ ] random(L), random vector inside a box of length L
// [ ] random_orient() random unit on unit sphere (marsaglia) 
#include<initializer_list> // inizializzazione vettori
#include<iostream> // input/output
#include<cmath>
#include<string>
// STRATEGIE PER RENDERE GENERICA LA CLASSE
//#define NT 3
// constexpr tells compiler that NT can be calculated at compile time
constexpr int NT=3;

//typedef float ntype; // as in C
using ntype=double; // use "using" to define an alias for the keyword "double"

class pvector
{
  ntype v[NT]; // private member
public:
  pvector() // void constructor
    {
      int i;
      for(i=0; i < NT; i++)
        {
          v[i]=0;
        }
    }

  // overloading of constructor for hangling curly brace initialization
  // ATTENTO le initializer list sono sempre dello stesso tipo (dichiarato dentro la classe con <> perche e di una libreria standard)
  pvector(std::initializer_list<ntype> list)  // overloading constructor
    {
      int c=0;
      // Auto NON e una lambda function, semplicemente inferisce sul tipo delle variabili. 
      // Invece che metterci double o ntype, ci ho messo qualcosa che lo fa in automatico
      for (auto el: list) // itero su tutti gli elementi della lista "list" e
                          // quindi "el" sarà
                          // ognuno di questi elementi 
        {
          if (c < NT)
            {
              v[c] = el;
            }
          c++;
        }
      for (;c < NT; c++) // gli elementi sono in numero di NT allora inizializzo 0
        {
          v[c]=0.0;
        }
    }


  ~pvector() // destructor
    {
    }

  // Note this method does not change data of object hence it is declared "const"
  void show(std::string s="") const // argomento di default 
    {
      std::cout << s << "(";
      for (int i=0; i < NT; i++)
        {
          std::cout << v[i];
          if (i < NT-1) 
            std::cout << ",";
        }
      std::cout << ")\n";
    }

  // since assign operator returns pvector objects by reference:
  // (A=B).set(0,2.3); it set element 0 o A to 2.3
  // Note this method changes data of object hence it is *not* declared "const"
  pvector& operator=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {

          (*this).v[i] = v2.v[i];
          // equivalently:
          // (*this).v[i] = v2.v[i];
        }
      return(*this);
    }



  // operator+: (*this)+v2
  // Note this method does not change data of object hence it is declared "const"
  pvector operator+(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = v[i] + v2.v[i];
          // equivalently:
          // vs.v[i] = (*this).v[i] + v2.v[i];
        }
      return vs;
    } 

  // operator-
  // Note this method does not change data of object hence it is declared "const"
  pvector operator-(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = v[i] - v2.v[i];
          // equivalently:
          // vs.v[i] = (*this).v[i] - v2.v[i];
        }
      return vs;
    } 


  // sum based on overloaded + operator (example of usage of (*this)) 
  // A.sum(B) questo è equivalente a scrivere A+B
  // Note this method does not change data of object hence it is declared "const"
  pvector sum(const pvector& v2) const
    {
      // *this represents the calling object (vector)
      // Questo è utilizzare l'ogetto stesso all'interno della classe, cioe, l'operatore + definito in questa classe? Boh non ho capito bene 
      return (*this)+v2;
    }


  // get
  // Note this method does not change data of object hence it is declared "const"
  ntype get(int i) const
    {
      return v[i];
      // equivalently:
      // return (*this).v[i];
    }

  // set
  ntype set(int i, ntype val) 
    {
      return v[i]=val;
      // equivalently: 
      // return (*this).v[i]=val;
    }

  // Esercizio da fare in laboratorio: implementare += e -= 
  pvector& operator+=(ntype num) 
    {
      for (int i=0; i < NT; i++)
        {
          (*this).v[i] += num;
          // equivalently:
          // vs.v[i] = (*this).v[i] - v2.v[i];
        }
      return (*this);
    } 
  pvector& operator-=(ntype num) 
    {
      for (int i=0; i < NT; i++)
        {
          (*this).v[i] -= num;
        }
      return (*this);
    } 

   // BTW le parentesi non sono obbligatorie
   pvector& operator+=(const pvector& v2) 
    {
      for (int i=0; i < NT; i++)
        {
          (*this).v[i] += v2.v[i];
          // equivalently:
          // vs.v[i] = (*this).v[i] - v2.v[i];
        }
      return (*this);
    } 
  pvector& operator-=(const pvector& v2) 
    {
      for (int i=0; i < NT; i++)
        {
          (*this).v[i] -= v2.v[i];
        }
      return (*this);
    } 
};
#endif