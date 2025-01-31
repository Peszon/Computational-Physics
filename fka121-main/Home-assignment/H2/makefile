# Compiler
CC = gcc

# Local installation directory
LOCAL=${HOME}/.local

# Compilation Flags
CFLAGS = \
     -pedantic \
     -Werror \
     -Wall \
     -g \
     -Iinclude \
     -Wl,-rpath,${LOCAL}/lib \
     -L${LOCAL}/lib \
     -I${LOCAL}/include \
     -fopenmp              # Added OpenMP flag for compilation

# Optimization Flags
CFLAGS_OPT = \
         -O2 \
         -march=native

# Libraries
LIB = \
      -lm \
      -lgsl \
      -lgslcblas \
      -fopenmp              # Added OpenMP flag for linking

# Object Files
OBJ = \
      obj/tools.o \
      obj/run.o \
      obj/lattice_utils.o  # Added lattice_utils object file

# Main Object File
MAIN = \
       obj/main.o

# Conditional Flags for Testing
ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

# Target for Building the Executable 'H2'
H2: obj _H2

# Linking the Object Files to Create 'H2'
_H2: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H2 $(LIB)

# Pattern Rule for Compiling Source Files to Object Files
obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

# Rule to Create 'obj' Directory if It Doesn't Exist
obj:
	mkdir -p obj

# Clean-Up Rules
clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f H2 run-test
	rm -rf obj obj_test

# Phony Targets
.PHONY: clean
