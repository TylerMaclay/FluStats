CC := g++
CCFLAGS := -Wall -Wextra -pedantic -std=c++17 -pthread

TARGETS := main
MAINS := $(addsuffix .o, $(TARGETS))
OBJ := FS_Stats.o OutputHandle.o config.o $(MAINS)
DEPS := OutputHandle.h FS_Stats.h config.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS)

$(TARGETS): % : $(filter-out $(MAINS), $(OBJ)) %.o
	$(CC) -o $@ $^ $(CCFLAGS)
