#include "tensors.hpp"
#include "config.hpp"
#include <iostream>
#include <initializer_list>
#include <string>

class Scalar {
    private: 
        ntype x;
    public: 
        Scalar() {
            x = 0;
        }
        
        Scalar(ntype num) {
            x = num;
        }

        ~Scalar() {}  
        // No overload is needed for this kind of variable 

        ntype get() const{ 
            return x; 
        }
};

class Vector {
    private: 
        ntype v[N1];
        int current_index;
    public:
        // Constructors 
        Vector() {
            int i;
            current_index = 0;
            for(i=0; i < N1; i++)
                {
                v[i]=0;
                }
        }

        Vector(std::initializer_list<ntype> list)  {
            int c=0;
            current_index = 0;
            for (auto el: list) 
                {
                if (c < N1)
                    {
                    v[c] = el;
                    }
                c++;
                }
            for (;c < N1; c++) // gli elementi sono in numero di NT allora inizializzo 0
                {
                v[c]=0.0;
                }
        }
        ~Vector() {}
        // Operators oveload and methods
        void show(std::string name = "") const {
            std::cout << name << "(";
            for (int i=0; i < N1; i++)
                {
                std::cout << v[i];
                if (i < N1-1) 
                    std::cout << ",";
                }
            std::cout << ")\n";
        }

        ntype get(int i) const {
            return (*this).v[i];
        }

        ntype set(int i, ntype val) {
            return (*this).v[i]=val;
        }

        double modulus() {
            return ((*this) * (*this));
        }
        // The imported vector is not modified so it's const 
        Vector& operator=(const Vector& v2)
        {
        for (int i=0; i < N1; i++)
            {
            (*this).v[i] = v2.v[i];
            }
        return(*this);
        }
        Vector operator+(const Vector& v2) {
            Vector v3;
            for (int i; i< N1; i++) {
                v3.v[i] = (*this).v[i] + v2.v[i]; 
            }
            return v3;
        }
        Vector operator-(const Vector& v2) {
            Vector v3;
            for (int i; i< N1; i++) {
                v3.v[i] = (*this).v[i] - v2.v[i]; 
            }
            return v3;
        }

        Vector& operator*(const Scalar& x) {
            for (int i; i< N1; i++) {
                this->v[i] = x.get() * this->v[i];
            }
            return (*this);
        }

        double operator*(const Vector& v2) {
            double prod = 0;
            for (int i; i< N1; i++) {
                prod += (*this).v[i] * v2.v[i];
            }
            return prod;
        }

        Vector operator^(const Vector& v2) const {
        if (N1 != 3) throw std::invalid_argument("X product is only defined for 3D vectors");
        Vector v3;
            v3.v[0] = v[1] * v2.v[2] - v[2] * v2.v[1];
            v3.v[1] = v[2] * v2.v[0] - v[0] * v2.v[2];
            v3.v[2] = v[0] * v2.v[1] - v[1] * v2.v[0];
        return v3 ;
    }

        Vector& operator+=(const Vector& v2) {
        for (int i=0; i < N1; i++)
            {
            (*this).v[i] += v2.v[i];
            }
        return (*this);
        } 

        Vector& operator-=(const Vector& v2) {
        for (int i=0; i < N1; i++)
            {
            (*this).v[i] -= v2.v[i];
            }
        return (*this);
        } 

        Vector& operator*=(const Scalar& x) {
        for (int i=0; i < N1; i++)
            {
            (*this).v[i] *= x.get();
            }
        return (*this);
        } 

        Vector& operator/=(const Scalar& x) {
        for (int i=0; i < N1; i++)
            {
            (*this).v[i] /= x.get();
            }
        return (*this);
        }  

        ntype& operator()(int index){
            return (*this).v[index];
        }

        friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "[";
        for (int i = 0; i < N1; ++i) {
            os << N1;
            if (i != N1 - 1) os << ", ";
        }
        os << "]";
        return os;
        }

        Vector& operator,(int val) {
            this->v[current_index++] = val;
            return *this;
        }

        bool operator==(const Vector& v2) {
            for (int i = 0; i < N1; ++i) {
                if (v2.v[i] != v2.v[i]) {
                    return false;
                }
            }
            return true;
        }
};