OBJS = main.o OneDNode.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -O -c
LFLAGS = -Wall -lm -O -o

OneDNode.o: OneDNode.cpp OneDNode.hpp
	$(CC) $(CFLAGS) OneDNode.cpp

main.o: main.cpp OneDNode.hpp OneDCell.hpp OneDGrid.hpp ExactFlux.hpp TimeMarching.hpp
	$(CC) $(CFLAGS) main.cpp

main: $(OBJS)
	$(CC) $(LFLAGS) main $(OBJS)

clean:
	\rm -f *.o *~ main