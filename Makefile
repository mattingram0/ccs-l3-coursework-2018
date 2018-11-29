#Constants that are subsituted into the file later when $() is used.
#CDFLAGS - compiler flags, LDFLAGS - library load flasg (-lm = load math), CC = compile command (gcc)
#Some useful constants:
#@ - Left hand side of the ':' in the recipe
#< - first item in the dependency list

CFLAGS = -O3 -march=native -D_GNU_SOURCE
PROFILE_FLAGS = -O3 -pthread -DLIKWID_PERFMON -I/ddn/apps/Cluster-Apps/likwid/4.1/include -L/ddn/apps/Cluster-Apps/likwid/4.1/lib
LIKWID_FLAGS = -llikwid -lm
GLIB_FLAGS = `pkg-config --cflags --libs glib-2.0`
LDFLAGS = -lm
CC = gcc

#Variable to test if we are profiling:
PROFILE = true

#List of object files, followed by the header files
OBJ = optimised-sparsemm.o basic-sparsemm.o utils.o
HEADER = utils.h

#Tells the make file not to do anything with files called clean, help, check, as these would conflict with our implicit rules
.PHONY: clean help check

#'make' with no arguments calls the first rule which is 'all' which calls 'sparemm'
all: sparsemm

#Help command
help:
	@echo "Available targets are"
	@echo "  clean: Remove all build artifacts"
	@echo "  check: Perform a simple test of your optimised routines"
	@echo "  sparsemm: Build the sparse matrix-matrix multiplication binary"

#Clean command, removes the executable and any object files
clean:
	-rm -f sparsemm $(OBJ)

check: sparsemm
	./sparsemm CHECK

#The sparsemm object file depends on sparsemm.c and all of the object files. If any of these files have changed, then the sparsemm executable will be remade
#How it is remade is defined by the recipe - so if any of the prerequisites have changed, and we run 'make sparsemm', we perform:
#gcc -03 -march=native -D_GNU_SOURCE -o sparsemm sparsemm.c optimised-sparsemm.o basic-sparsemm.o utils.o
#which recompiles sparsemm.c optimised.o, basic-sparsemm.o and utils.o into our sparsemm executable
sparsemm: sparsemm.c $(OBJ)
	$(CC) $(CFLAGS) $(GLIB_FLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

profile: sparsemm.c $(OBJ)
	$(CC) $(PROFILE_FLAGS) $(GLIB_FLAGS) -o $@ $< $(OBJ) $(LIKWID_FLAGS)

#As our make file works recursively, if when we run 'make sparsemm' our other .c files (basic-sparsemm.c etc) have changed, then our .o should be remade too. This recipe handles this:
#Our target is our %.o file, and its prerequisites/dependencies are our %.c files and the header files. If a .c (or the header file) has changed, then we perform:
#gcc -O3 -march=native  -D_GNU_SOURCE -c -o %.o %.c
%.o: %.c $(HEADER)
ifeq ($(PROFILE),true)
	$(CC) $(PROFILE_FLAGS) $(GLIB_FLAGS) $< -o $@ -c
else
	$(CC) $(CFLAGS) $(GLIB_FLAGS) $< -o $@ -c
endif
