#===General Variables===#
CC=gcc
#CFLAGS=-Wall -Wextra -g3 -Ofast
CFLAGS=-Wall -Wextra -g3

all: makeAll

makeAll: makeHelper makeRandom makeStats  makeAlgebra makeMachine makeMain
	$(CC) $(CFLAGS) helper.o random.o stats.o  algebra.o machine_learning.o main.o -o machine_learning -ldl -lm -lblas -llapack -lpthread 

makeMain: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o 
	#$(CC) $(CFLAGS) -c `pkg-config --cflags opencv` main.c -o main.o 

makeMachine: machine_learning.c machine_learning.h
	$(CC) $(CFLAGS) -c machine_learning.c -o machine_learning.o 

makeAlgebra: linear_algebra.c linear_algebra.h
	$(CC) $(CFLAGS) -c linear_algebra.c -o algebra.o 

makeStats: stats.c stats.h
	$(CC) $(CFLAGS) -c stats.c -o stats.o 

makeRandom: random_numbers.c random_numbers.h
	$(CC) $(CFLAGS) -c random_numbers.c -o random.o 

makeHelper: helper.c helper.h macros.h
	$(CC) $(CFLAGS) -c helper.c -o helper.o


.PHONY: clean

clean:
	rm -f *~ $(ODIR)/*.o $(IDIR)/*~
