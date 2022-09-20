PARA = -std=c++14 -Wall -O3 -g -DNDEBUG
com: graph.cpp main.cpp
	g++ -c graph.cpp -o graph.o $(PARA)	
	g++ main.cpp graph.o -o fasthare  $(PARA)
