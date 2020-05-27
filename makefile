CFLAGS = -std=c++11 -Wall -O3 -c -g -DDMEASURE #-DPERTURB_CRIT #-DDEBUG #-DDEBUG_CHI_COMPUTATION

.PHONY = all

all: subgraph clean

subgraph: main.o query_graph.o input_graph.o vertex.o node.o
	g++ -o subgraph -Wall -O3 main.o query_graph.o input_graph.o vertex.o node.o

node.o: node.cpp node.h
	g++ $(CFLAGS) node.cpp

vertex.o: vertex.cpp vertex.h node.h const.h query_graph.h
	g++ $(CFLAGS) vertex.cpp

query_graph.o: query_graph.cpp query_graph.h node.h
	g++ $(CFLAGS) query_graph.cpp

input_graph.o: input_graph.cpp input_graph.h vertex.h node.h const.h
	g++ $(CFLAGS) input_graph.cpp

main.o: main.cpp query_graph.h input_graph.h vertex.h node.h const.h
	g++ $(CFLAGS) main.cpp


clean:
	rm *.o
