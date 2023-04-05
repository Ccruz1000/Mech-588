run: src/Roe.o src/Roe.h src/main.o src/test.o src/Steger.o src/functions.o
	g++ -g src/Steger.o src/functions.o src/test.o src/Roe.o -o test
	./test input.txt
	g++ -g src/Steger.o src/functions.o src/Roe.o src/main.o -o run

all: src/Roe.cc src/Roe.h src/main.cc src/test.cc src/Steger.o src/functions.o
	g++ -g -c src/Roe.cc -o src/Roe.o
	g++ -g -c src/Steger.cc -o src/Steger.o
	g++ -g -c src/functions.cc -o src/functions.o
	g++ -g -c src/main.cc -o src/main.o
	g++ -g -c src/test.cc -o src/test.o
	g++ -g src/Steger.o src/functions.o src/test.o src/Roe.o -o test
	./test input.txt
	g++ -g src/Roe.o src/Steger.o src/functions.o src/main.o -o run

test: src/test.o src/Roe.o
	g++ -g src/test.o src/Roe.o src/Steger.o src/functions.o -o test

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

src/test.o: src/test.cc
	g++ -g -c src/test.cc -o src/test.o

src/Roe.o: src/Roe.cc src/Roe.h
	g++ -g -c src/Roe.cc -o src/Roe.o

src/Steger.o: src/Steger.cc src/Steger.h
	g++ -g -c src/Steger.cc -o src/Steger.o

src/functions.o: src/functions.cc src/Steger.h
	g++ -g -c src/functions.cc -o src/functions.o

clean:
	rm src/*.o run test