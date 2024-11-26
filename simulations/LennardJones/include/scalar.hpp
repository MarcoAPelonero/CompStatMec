#ifndef SCALAR_HPP
#define SCALAR_HPP

#include "config.hpp"

class Scalar {
private:
    ntype x;
public:
    Scalar();
    Scalar(ntype num);
    ~Scalar(); 

    ntype get() const;
};

#endif 