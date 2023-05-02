run: src/poisson.o  src/poisson.h src/main.o src/poisson4th.h src/poisson4th.o
	g++ -g src/poisson.o src/main.o src/poisson4th.o -o run -lm

src/poisson.o: src/poisson.cc
	g++ -g -c src/poisson.cc -o src/poisson.o

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

src/poisson4th.o: src/poisson4th.cc
	g++ -g -c src/poisson4th.cc -o src/poisson4th.o

clean:
	rm src/*.o run
