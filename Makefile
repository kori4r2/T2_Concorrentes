CFLAGS = -g -Wall -std=c++11
LINKFLAGS = -pthread -std=c++11
PROJECT = GaussJordan
CC = g++
DEBUGDIR = tests
BINDIR = bin
SRCDIR = src
LIBDIR = lib
LIBS := $(wildcard $(LIBDIR)/*hpp)
SRCS := $(wildcard $(SRCDIR)/*cpp)
OBJS := $(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o, $(SRCS))

#TO DO: Run on remote server and test it

all : build

build : $(BINDIR)
	$(CC) $(LINKFLAGS) -o $(SERVER) $(SERVEROBJS)

$(DEBUGDIR) :
	mkdir -p $(DEBUGDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

$(BINDIR)/%.o : $(SRCDIR)/%.cpp $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

hostfile :
	scp -P 22200 hosts usuario@halley.lasdpc.icmc.usp.br:/home/usuario/ #TO DO

clean :
	rm -rf $(BINDIR)
	rm -rf $(DEBUGDIR)
	rm -f $(PROJECT).zip
	rm -f $(PROJECT)
	rm -f debug*.txt
	clear

run :
	./$(PROJECT)

.zip : clean
	zip $(PROJECT).zip $(SRCS) $(LIBS) Makefile *.pdf Authors.txt

#debugServer: $(DEBUGDIR) all
#	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(PROJECT) > $(DEBUGDIR)/output.txt 2> $(DEBUGDIR)/error.txt
