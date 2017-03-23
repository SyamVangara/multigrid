#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
CC = gcc
CFLAGS  = -Wall -std=c99

# External libraries:
LIBS = -lm

# Pre-defined macros for conditional compilation
DEFS = -DDEBUG_FLAG -DEXPERIMENTAL=0

BIN = poisson

DEP = header.h

OBJS = solver.o mesh.o array.o

$(BIN): $(OBJS) $(DEP)
	$(CC) $(CFLAGS) $(DEFS) $(OBJS) $(BIN).c -o $(BIN) $(LIBS)

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $(DEFS) $< -o $@

clean: 
	$(RM) count *.o *~

depend:
	makedepend -Y -- $(CFLAGS) $(DEFS) -- *.c
