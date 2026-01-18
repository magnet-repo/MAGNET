CC=gcc
CFLAGS=-O3 -march=native -Wall -I.
LDFLAGS = -lm

OBJS=cpucycles.o fastrandombytes.o gaussian_avx.o hash.o keygen.o mint.o poly_mult.o poly_q.o randombytes.o sammat.o sampleb.o setup.o spend.o test.o verify.o utils.o params.o gadgets.o masked_algorithms.o base_gadgets.o sign_gadgets.o params-GR19.o utils-GR19.o

PROGRAM=ringct

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROGRAM) $(LDFLAGS) libXKCP.a

clean:;
	rm -f *.o $(PROGRAM)
