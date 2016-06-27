CC=g++
CFLAGS=-g -Wall -O2 -fopenmp 
LDFLAGS=-g -Wall -c -fopenmp
EXECUTABLES=calc
SOURCES=coef.cpp info.cpp pref_co.cpp file_pro.cpp spatial.cpp concentration.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o $(EXECUTABLES)

.cpp.o:
	$(CC) $(LDFLAGS) $(SOURCES)

clean:
	rm *.o 

