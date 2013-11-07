CC = g++
CFLAGS = -Wall -O3 -fstrict-aliasing -std=c++0x
CFLAGS += $(shell pkg-config --cflags eigen3)
LIBS += $(shell pkg-config --libs eigen3)
LDFLAGS =
OBJECTS = lattice.o driver.o
PROGRAM = lbm.out

COMPILE = $(CC) $(CFLAGS)
LINK 		= $(CC) $(LDFALGS)

$(PROGRAM): $(OBJECTS)
	$(LINK) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o:%.cpp
	$(COMPILE) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.exe *.o *.out 
