CC = g++
CFLAGS = -Wall -std=c++11 -g -O3
CFLAGS += $(shell pkg-config --cflags eigen3)
LIBS += $(shell pkg-config --libs eigen3)
LDFLAGS =
OBJECTS = lattice.o node.o driver.o
PROGRAM = lbm.out

COMPILE = $(CC) $(CFLAGS) $(LDFLAGS)
LINK 		= $(CC) $(CFLAGS) $(LDFLAGS)

$(PROGRAM): $(OBJECTS)
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o:%.cpp
	$(COMPILE) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.exe *.o *.out 
