CC = gcc
CXX = g++

CFLAGS = -g -W -Wall -O0 -m32
LDFLAGS = -m32

HEADERS =
SOURCES = setup.cc sign.cc verify.cc evaluate.cc
LIBRARIES =
OBJECTS =  
TARGETS = setup sign verify evaluate
INCNTL = -I/ #NTL include dir
LIBNTL = -L/ #NTL lib dir
INCBOTAN = -I/ #Botan include dir
LIBBOTAN = -L/ #Botan lib dir

default: $(TARGETS)

setup : setup.cc
	$(CXX) $(CFLAGS) -o $@ $^ $(INCNTL) $(LIBNTL) $(LFLAGS) -lntl -lm
sign : sign.cc
	$(CXX) $(CFLAGS) -o $@ $^ $(INCBOTAN) $(LIBBOTAN) $(INCNTL) $(LIBNTL) $(LFLAGS) -lntl -lbotan-1.10 -lm

verify : verify.cc
	$(CXX) $(CFLAGS) -o $@ $^ $(INCBOTAN) $(LIBBOTAN) $(INCNTL) $(LIBNTL) $(LFLAGS) -lntl -lbotan-1.10 -lm

evaluate : evaluate.cc
	$(CXX) $(CFLAGS) -o $@ $^ $(INCBOTAN) $(LIBBOTAN) $(LIBNTL) $(INCNTL) $(LFLAGS) -lntl -lbotan-1.10 -lm

Makefile.dependencies:: $(SOURCES) $(HEADERS) 
	$(CC) $(CFLAGS) -MM $(SOURCES) > Makefile.dependencies
-include Makefile.dependencies

clean:
	@rm -f $(TARGETS) *.o 

