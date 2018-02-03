
CC=g++
LD=g++ -fopenmp
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g -std=c++14 -O3 

OBJS = main.o scanner.o swz_len.o swz_3rrs.o 
DEPS = scanner.h swz_len.h swz_3rrs.h

MAIN = out

$(MAIN):$(OBJS)
	$(LD) -pg $(OBJS) $(CFLAGS) -o $(MAIN) 

%.o: %.cpp $(DEPS)
	$(LD) $(CFLAGS) -c $< -pg

clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res
	$(RM)  *.o