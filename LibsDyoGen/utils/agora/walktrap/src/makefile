CC = g++
PYTHON_INCLUDES = $(shell python-config --includes)
CFLAGS = -W -Wall -O3 -static -fPIC
LDFLAGS = -shared -fPIC

all : walktrap _walktrap.so

walktrap : walktrap.o communities.o graph.o heap.o
	$(CC) -o $@ $^

_walktrap.so : pywalktrap.o communities.o graph.o heap.o
	$(CC) -o $@ $^ $(LDFLAGS)

py%.o : py%.cpp
	$(CC) -c $< $(CFLAGS) $(PYTHON_INCLUDES)

%.o : %.cpp
	$(CC) -c $< $(CFLAGS)
clean :
	rm -f *.o *.pyo

moreclean :
	rm -f *.o *.pyo _walktrap.so walktrap

