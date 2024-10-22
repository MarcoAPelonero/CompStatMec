#include <iostream>
#include <vector>

class Rectangle {
    // All'interno di private in genere ci vanno solo le variabili di classe 
    // e i metodi che servono solamente alla classe stessa
    private:
        double height, width;
    public:
        // Tra i metodi public c'e sempre il constructor, un metodo chiamato
        //ESATTAMENTE come la classe, che è il metodo __init__ di python
        Rectangle(double h, double w) : height(h), width(w) {}
        // In qeusto caso il constructor non deve effettivmente processare niente ma solo assegnare le variabili 
        //ai double giusti nel private, quindi mettiamo le parentesi e basta
        /*
        double area() {
            double area = height * width;
            return area;
        }
        */
       // Dichiarazione piu efficiente del metodo area, che restituisce l'area del rettangolo
       double area() {
            return height * width;
        }
};

// Esempio di oveloading degli overatori

class Vector {
    private: 
        int x, y;
    public:
        Vector (int x, int y) : x(x), y(y) {}

        // Operator è una reoutine che sostanzialmente prende l'operatore
        // che gli metti davanti, e quando dell'operatore viene usato con un elemento della classe,
        // allora viene sovrascritto da questa cosa.
        Vector operator + (const Vector& v) {
            return Vector(x + v.x, y + v.y);
        }
        // Since << is not part of Vector class, we add firend, so it may access private variables 
        // Facciamo l'overload di un operatore di out cosi possiam osta,mparae il vettore con estrema facilità

        friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
            os << "Vector: (" << v.x << ", " << v.y << ")" << std::endl;
            return os;
        }
        void display() {
            std::cout << "Vector: (" << x << ", " << y << ")" << std::endl;
        }
        
};// Adesso vediamo la composizione di classi 

class Engine {
    public:
    void start() {
        std::cout << "Engine Started" << std::endl;
    }
};

class Car {
    private:
        Engine engine;
    public:
        void start() {
            engine.start();
            std::cout << "Car started" << std::endl;
        }
};

// Adesso vediamo come le classi ereditano e come
// si possono overloadare anche i loro metodi dopo l'eredità

class Shape {
    // So that the classes that may inherit Shape can use these variables
    protected: 
        double width, height;
    public: 
        Shape(double width, double height) : width(width), height(height) {}

        // I plan to override this in future classes, so i make it a virtual method
        virtual double area() {
            return width * height;
        }
};
// Il metterci la classe Shape accanto segna l'ereditarieta
class Rect2: public Shape {
    public: 
    Rect2(double w, double h) : Shape(w,h) {}

    double area() override {
        return width * height;
    }
};
// Vediamo come ereditare Shape, e qualcosa di utile anche per altre figure geometriche
class Triangle : public Shape {
public:
    // Constructor calls the base class constructor
    Triangle(double w, double h) : Shape(w, h) {}

    // Overriding the area method for Triangle
    double area() override {
        return (width * height) / 2;
    }
};

// Per esempio posso fare invece un function overload in certi casi e
//d efinire i casi di uso per specifiche funzioni 

class Print {
    public:
     void display(int value) {
        std::cout << "Integer: " << value << std::endl;
    }

    // Overload for double
    void display(double value) {
        std::cout << "Double: " << value << std::endl;
    }

    // Overload for string
    void display(const std::string& value) {
        std::cout << "String: " << value << std::endl;
    }
};

// usando le classi e pratica comune utilizzare delle funzioni virtuali
// in modo tale da poter ricostruire il path di arrivo di determinate funzioni e metodi

class Animal {
    public: 
    virtual void sound() {
        std::cout << "The animal makes a sound:" << std::endl;
    }
};

class Dog: public Animal {
    public: 
        void sound() override {
            std::cout << "Dog barks" << std::endl;
        }
};

class Cat: public Animal {
    public: 
        void sound() override {
            std::cout << "Cat meows" << std::endl;
        }
};