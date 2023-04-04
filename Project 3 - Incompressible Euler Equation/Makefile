run: src/Roe.o src/Roe.h src/main.o src/test.o
	g++ -g src/test.o src/Roe.o src/Roe.h -o test
	./test input.txt
	g++ -g src/Roe.o src/main.o src/Roe.h -o run

all: src/Roe.cc src/Roe.h src/main.cc src/test.cc
	g++ -g -c src/Roe.cc -o src/Roe.o
	g++ -g -c src/main.cc -o src/main.o
	g++ -g -c src/test.cc -o src/test.o
	g++ -g src/test.o src/Roe.o src/Roe.h -o test
	./test input.txt
	g++ -g src/Roe.o src/Roe.h src/main.o -o run

test: src/test.o src/Roe.o
	g++ -g src/test.o src/Roe.o src/Roe.h -o test

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

src/test.o: src/test.cc
	g++ -g -c src/test.cc -o src/test.o

src/Roe.o: src/Roe.cc
	g++ -g -c src/Roe.cc -o src/Roe.o

clean:
	rm src/*.o run test