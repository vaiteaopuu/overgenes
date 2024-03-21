CC=gcc
BIN=./bin

OBJS2 = \
	global.o\
	seq_utils.o

OBJS1 = \
	design.o\
	seq_utils.o

OBJS3 = \
	mcmc.o\
	seq_utils.o

OBJS4 = \
	local.o\
	seq_utils.o

CFLAGS+= -g
# CFLAGS+=  -O2 -flto -Wall -Wextra -DHAVE_INLINE
# CFLAGS+= -Ofast -march=native -flto -Wall -Wextra -DHAVE_INLINE

all: design global mcmc local

design: $(OBJS1)
	$(CC) $(CFLAGS) -o $(BIN)/design $(OBJS1) -lm

global: $(OBJS2)
	$(CC) $(CFLAGS)  -o $(BIN)/global $(OBJS2) -lm

mcmc: $(OBJS3)
	$(CC) $(CFLAGS)  -o $(BIN)/mcmc $(OBJS3) -lm

local: $(OBJS4)
	$(CC) $(CFLAGS)  -o $(BIN)/local $(OBJS4) -lm

clean:
	rm -rf *.o $(BIN)/global $(BIN)/design  $(BIN)/local $(BIN)/mcmc
