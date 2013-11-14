CC = g++
CFLAGS = -Wall -O3 -std=c++11 -g
CFLAGS += $(shell pkg-config --cflags eigen3)
LIBS += $(shell pkg-config --libs eigen3)
LDFLAGS =
OBJECTS = lattice.o driver.o
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
