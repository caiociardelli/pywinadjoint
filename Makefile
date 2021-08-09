# Makefile for PyWinEPAdjoint

# Compiler and linker
CC := gcc
# Flags for compiler 1
FLAGS1 := \
	-std=c99 \
	-Wall \
	-Wextra \
	-shared \
	-fPIC \
	-O3
# Flags for compiler 2
FLAGS2 := \
	-std=c99 \
	-Wall \
	-Wextra \
	-O3
# Source directory
SRC := src
# Includes
INC := -I$(SRC)/headers
# Objects directory
OBJ := obj
# Libraries 1
LIB1 := -lfftw3 -lm
# Libraries 2
LIB2 := -lm
# Binaries directory 1
BIN1 := lib
# Binaries directory 2
BIN2 := bin

# Default targets 1
DEFAULT1 := \
	libsignal.so \
	libdistance.so \
	libweights.so

# Default targets 2
DEFAULT2 := expand

# Complete path for binaries 1
DFT1 := $(patsubst %, $(BIN1)/%, $(DEFAULT1))

# Complete path for binaries 2
DFT2 := $(patsubst %, $(BIN2)/%, $(DEFAULT2))

# Command used for cleaning
RM := rm -rf
 
#
# Compilation and linking
#
all: objDirectory binDirectory $(DFT1) $(DFT2)
	@ echo 'Finished building binary!'

$(BIN1)/libsignal.so: $(OBJ)/signal.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS1) $(INC) $^ -o $@ $(LIB1)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN1)/libdistance.so: $(OBJ)/distance.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS1) $(INC) $^ -o $@ $(LIB2)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN1)/libweights.so: $(OBJ)/weights.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS1) $(INC) $^ -o $@ $(LIB2)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN2)/expand: $(OBJ)/expand.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS2) $^ -o $@
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(OBJ)/signal.o: $(SRC)/signal.c
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS1) $(INC) -c $(SRC)/signal.c -o $@
	@ echo ' '

$(OBJ)/distance.o: $(SRC)/distance.c
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS1) $(INC) -c $(SRC)/distance.c -o $@
	@ echo ' '

$(OBJ)/weights.o: $(SRC)/weights.c
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS1) $(INC) -c $(SRC)/weights.c -o $@
	@ echo ' '

$(OBJ)/expand.o: $(SRC)/expand.c
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS2) -c $(SRC)/expand.c -o $@
	@ echo 'Finished building binary: $@'
	@ echo ' '

objDirectory:
	@ mkdir -p $(OBJ)

binDirectory:
	@ mkdir -p $(BIN1) $(BIN2)

clean:
	$(RM) $(OBJ)/ $(BIN1)/ $(BIN2)/

.PHONY: all clean
