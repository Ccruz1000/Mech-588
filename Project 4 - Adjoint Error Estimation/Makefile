run: src/poisson.o  src/poisson.h src/main.o
	g++ -g src/poisson.o src/main.o -o run -lm

src/poisson.o: src/poisson.cc
	g++ -g -c src/poisson.cc -o src/poisson.o

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

clean:
	rm src/*.o run
