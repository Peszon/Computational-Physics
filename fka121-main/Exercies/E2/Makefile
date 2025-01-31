CC ?= gcc # c-compiler
CFLAGS = # compiler flags such as -O3 -W
LIBS = -lgsl -lgslcblas -lm # libraries such as -lm -lgsl -lgslcblas

PROGRAM = E2 # name of executable to be built

LOCAL = ${HOME}/.local
CFLAGS += \
	-Wl,-rpath,${LOCAL}/lib \
	-I${LOCAL}/include \
	-I/home/felixpersson/kompfys/fka121-main/Exercies/E2/include \
	-L${LOCAL}/lib \
	-Werror \
	-pedantic \
	-g \
	-O0 \
	-fsanitize=address \
	-fno-omit-frame-pointer

# Specify the source files
SOURCES = src/main.c src/tools.c
OBJECTS = $(SOURCES:.c=.o)

# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

# compile program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

# clean up object files and executable
clean:
	rm -f $(OBJECTS) $(PROGRAM)