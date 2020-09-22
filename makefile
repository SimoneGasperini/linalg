all: linalg.a clean

linalg.a: Polynomial.o Vector.o Matrix.o linalg.o
	ar rvs linalg.a Polynomial.o Vector.o Matrix.o linalg.o

linalg.o: linalg.cpp
	g++ -c linalg.cpp

Matrix.o: Matrix.cpp
	g++ -c Matrix.cpp

Vector.o: Vector.cpp
	g++ -c Vector.cpp

Polynomial.o: Polynomial.cpp
	g++ -c Polynomial.cpp

clean:
	rm -f *.o