run: src/mesh.o src/main.o src/field.o src/test.o
	g++ -g src/test.o src/mesh.o -o test
	./test
	g++ -g src/mesh.o src/main.o src/field.o -o run

all: src/Mesh_Read.cc src/main.cc src/field_functions.cc
	g++ -g -c src/Mesh_Read.cc -o src/mesh.o
	g++ -g -c src/main.cc -o src/main.o
	g++ -g -c src/field_functions.cc -o src/field.o
	g++ -g -c src/test.cc -o src/test.o
	g++ -g src/test.o src/mesh.o -o test
	./test
	g++ -g src/mesh.o src/main.o src/field.o -o run

test: src/test.o
	g++ -g src/test.o -o test

src/mesh.o: src/Mesh_Read.cc
	g++ -g -c src/Mesh_Read.cc -o src/mesh.o

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

src/field.o: src/field_functions.cc
	g++ -g -c src/field_functions.cc -o src/field.o

src/test.o: src/test.cc
	g++ -g -c src/test.cc -o src/test.o

clean:
	rm src/*.o run test
