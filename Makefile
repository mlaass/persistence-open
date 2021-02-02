all:
	g++ -I./include/  -std=c++11 -ofast -march=native  -o main main.cpp -g
fast:
	g++ -I./include/  -std=c++11 -o4 -march=native  -o main main.cpp
