g++ -c sde_mult.cpp
#g++ -g -O0 sde_mult.cpp -lgsl -lgslcblas -lm -lstdc++ -ltrapfpe
g++  sde_mult.o -lgsl -lgslcblas -lm -lstdc++ -ltrapfpe
