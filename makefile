all: mst-backbone test

mst-backbone: mst-backbone.o
	g++ -o mst-backbone mst-backbone.o -lstdc++fs

all: mst-backbone

mst-backbone.o: mst-backbone.cpp
	g++ -c mst-backbone.cpp -o mst-backbone.o -O3 -I. -std=c++11

clean:
	rm mst-backbone mst-backbone.o


