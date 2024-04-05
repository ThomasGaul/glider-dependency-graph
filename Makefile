main: main.o random.o LtL.o Geometry.o DependencyGraph.o
	g++ -o main main.o random.o LtL.o Geometry.o DependencyGraph.o

main.o: main.cpp LtL.h DependencyGraph.h Geometry.h
	g++ -c -O3 main.cpp


random.o: random.cpp random.h
	g++ -c -O3 random.cpp

LtL.o: LtL.cpp LtL.h random.h
	g++ -c -O3 LtL.cpp

Geometry.o: Geometry.cpp Geometry.h
	g++ -c -O3 Geometry.cpp

DependencyGraph.o: DependencyGraph.cpp DependencyGraph.h
	g++ -c -O3 DependencyGraph.cpp

clean:
	rm *.o main