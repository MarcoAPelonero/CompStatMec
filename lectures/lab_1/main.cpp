#include <iostream>
#include <vector>
#include "classes.hpp"

// Si puo anche aggiungere un namespace, all'interno del quale ci sono
// delle operazioni custom, racchiuse all'interno del namespace
// questo diventa utile quando i programmi diventano enormi

namespace MathOperations {
    int add(int a, int b) {
        return a + b;
    }

    int multiply(int a, int b) {
        return a * b;
    }
}

int main() {
    Rectangle rect(5.0, 3.0);
    std::cout << "Rectangle area: " << rect.area() << std::endl;

    Vector v1(1,2);
    Vector v2(2,4);
    // Come vedi il + adesso è illuminato ed e perche verra compilato come
    // v1.operator+(v2)
    Vector v3 = v1 + v2;
    v3.display();
    // Posso anche usare l'overloading dell'op <<
    std::cout << v1 << v2 << v3 << std::endl;

    Car mazda;
    mazda.start();

    std::cout << "Sum: " << MathOperations::add(5, 3) << std::endl;  // Output: Sum: 8
    std::cout << "Product: " << MathOperations::multiply(5, 3) << std::endl;  // Output: Product: 15

    Rect2 rect2(10, 5);
    Triangle tri(10, 5);

    std::cout << "Area of Rectangle: " << rect2.area() << std::endl;  // Output: 50
    std::cout << "Area of Triangle: " << tri.area() << std::endl;

    Print printer;
    
    printer.display(5);             // Calls the int version, Output: Integer: 5
    printer.display(3.14);          // Calls the double version, Output: Double: 3.14
    printer.display("Hello");

    // Creo un puntatore a classe Animal
    // e dopo sovrascrivo i suoi metodi con quelli di dog e cat.
    // In questi casi per utilizzare gli attriuti non si usa il . ma il ->
    Animal* animalPtr;

    Dog dog;
    Cat cat;

    animalPtr = &dog;
    animalPtr->sound();  // Output: Dog barks (calls Dog's sound method)

    animalPtr = &cat;
    animalPtr->sound();
}