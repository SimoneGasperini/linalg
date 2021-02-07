cmake_minimum_required(VERSION 3.16)

project(linalg)


add_executable(

    example.out

    example.cpp

    linalg.cpp
    Vector.cpp
    Polynomial.cpp
    Matrix.cpp

)