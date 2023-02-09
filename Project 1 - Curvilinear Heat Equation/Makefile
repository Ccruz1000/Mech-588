run: main.o Fcurv.o
	g++ -g main.o Fcurv.o -o run

test: test.o Fcurv.o
	g++ -g test.o Fcurv.o -o test

main.o: main.cc
	g++ -g -c main.cc

Fcurv.o: Fcurv.cc
	g++ -g -c Fcurv.cc

test.o: test.cc
	g++ -g -c test.cc

clean:
	rm *.o test run

