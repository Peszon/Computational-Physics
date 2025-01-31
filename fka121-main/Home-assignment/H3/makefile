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
     -fopenmp \
     -fsanitize=address \
     -fno-omit-frame-pointer            

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
      obj/H3_utils.o

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

# Target for Building the Executable 'H3'
H3: obj _H3

# Linking the Object Files to Create 'H3'
_H3: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H3 $(LIB)

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
	rm -f H3 run-test
	rm -rf obj obj_test

# Phony Targets
.PHONY: clean
