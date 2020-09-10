all: linalg.a clean

linalg.a: Polynomial.o Matrix.o Vector.o linalg.o
	ar rvs linalg.a Polynomial.o Matrix.o Vector.o linalg.o

linalg.o: linalg.cpp
	g++ -c linalg.cpp

Vector.o: Vector.cpp
	g++ -c Vector.cpp

Matrix.o: Matrix.cpp
	g++ -c Matrix.cpp

Polynomial.o: Polynomial.cpp
	g++ -c Polynomial.cpp

clean:
	rm -f *.o