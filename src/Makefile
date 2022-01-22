all: tumor

tumor: objects
	g++ -o tumor main.o grid.o cell.o -lglut -lGLU -lGL

objects:
	g++ -c grid.cpp -o grid.o
	g++ -c cell.cpp -o cell.o
	g++ -c main.cpp -o main.o

clean: 
	rm *.o tumor