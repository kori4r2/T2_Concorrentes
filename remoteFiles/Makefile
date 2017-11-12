CFLAGS = -g -O3 -fopenmp
LINKFLAGS = -fopenmp -lm
PROJECT = GaussJordan
CC = mpicc
BINDIR = bin
SRCDIR = src
LIBDIR = lib
LIBS := $(wildcard $(LIBDIR)/*h)
SRCS := $(wildcard $(SRCDIR)/*c)
OBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SRCS))
all : unzip build

unzip :
	unzip GaussJordan.zip

run: build
	mpirun --hostfile hosts -np 12 GaussJordan > resultado.txt

build : $(BINDIR) $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(PROJECT)

$(BINDIR) :
	mkdir -p $(BINDIR)

$(BINDIR)/%.o : $(SRCDIR)/%.c $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

clean :
	rm -rf $(BINDIR)
	rm -rf $(DEBUGDIR)
	rm -f $(PROJECT).zip
	rm -f $(PROJECT)
	rm -f debug*.txt
	clear

