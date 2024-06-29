CC = mpic++
CFLAGS = -Wall -O2
LDFLAGS = -lmpi

SRCS = poisson.cpp methods.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = t1

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@ -I.

clean:
	rm -f $(OBJS) $(EXEC)

run:
	mpiexec -n $(shell read -p "Enter the number of processes: " nproc; echo $$nproc) ./$(EXEC)
