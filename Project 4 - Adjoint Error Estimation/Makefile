run: src/poisson.o  src/poisson.h src/main.o src/poisson4th.h src/poisson4th.o src/Q_Ans.o src/Q_Ans.h
	g++ -g src/poisson.o src/main.o src/poisson4th.o src/Q_Ans.o -o run -lm

src/poisson.o: src/poisson.cc
	g++ -g -c src/poisson.cc -o src/poisson.o

src/main.o: src/main.cc
	g++ -g -c src/main.cc -o src/main.o

src/poisson4th.o: src/poisson4th.cc
	g++ -g -c src/poisson4th.cc -o src/poisson4th.o

src/Q_Ans.o: src/Q_Ans.cc
	g++ -g -c src/Q_Ans.cc -o src/Q_Ans.o

clean:
	rm src/*.o run
