CFLAGS = -m64 -Ofast -flto -funroll-loops -Wall -Wextra -ansi -pedantic -I.. -I/usr/local/include/dnest3
LIBS =  -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_system -lboost_thread -L..

OBJECTS = Data.o main.o MyDistribution.o MyModel.o radial_v.o

main: $(OBJECTS)
	g++ -o $@ *.o $(LIBS)
	# rm -f *.o

%.o: %.cpp
	g++ $(CFLAGS) -c $^

%.o: %.f90
	gfortran -c $^

clean:
	rm main

cleanall:
	rm -f *.o