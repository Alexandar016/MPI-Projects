# Compiler
FC = mpif90

# Compiler flags

CFLAG = -O2 -g

# Source files
SRCS = jacobi.f90
OBJS = $(SRCS:.f90=.o)

# Executable name
EXEC = jacobi

# Default target
all: $(EXEC)

# Rule to link the executable
$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS)

# Rule to compile Fortran source files
%.o: %.f90
	$(FC) -c $< -o $@

# Run target
run:
	@read -p "Enter number of processes (x) to run: " PROCESSES; \
	mpirun -np $$PROCESSES ./$(EXEC)

# Clean target
clean:
	rm -f $(EXEC) $(OBJS)