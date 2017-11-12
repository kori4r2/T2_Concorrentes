CFLAGS = -g -Wall -O3 -fopenmp
LINKFLAGS = -fopenmp -lm
PROJECT = GaussJordan
CC = mpicc
DEBUGDIR = tests
BINDIR = bin
SRCDIR = src
LIBDIR = lib
SERVERFILES := $(wildcard remoteFiles/*)
LIBS := $(wildcard $(LIBDIR)/*h)
SRCS := $(wildcard $(SRCDIR)/*c)
OBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SRCS))

all : build

build : $(BINDIR) $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(PROJECT)

$(DEBUGDIR) :
	mkdir -p $(DEBUGDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

$(BINDIR)/%.o : $(SRCDIR)/%.c $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

remote : sendfiles
	ssh gpra07@halley.lasdpc.icmc.usp.br -p 22200

sendfiles : .zip
	scp -P 22200 $(PROJECT).zip $(SERVERFILES) gpra07@halley.lasdpc.icmc.usp.br:/home/gpra07/

clean :
	rm -rf $(BINDIR)
	rm -rf $(DEBUGDIR)
	rm -f $(PROJECT).zip
	rm -f $(PROJECT)
	rm -f debug*.txt
	clear

run : build
	./$(PROJECT)

mpi_run: build
	mpirun -np 2 $(PROJECT)

.zip : clean
	zip $(PROJECT).zip $(SRCS) $(LIBS) Makefile

debug: $(DEBUGDIR) all
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes mpirun -np 2 $(PROJECT) > $(DEBUGDIR)/output.txt 2> $(DEBUGDIR)/error.txt
	diff resultadoCorreto.txt $(DEBUGDIR)/output.txt > $(DEBUGDIR)/diff.txt
