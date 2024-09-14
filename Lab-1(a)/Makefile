# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -Wall -g

# Linker flags for LAPACK and BLAS (Update with correct path)
LDFLAGS = -L/opt/homebrew/Cellar/openblas/0.3.28/lib -lopenblas

# Target executable
TARGET = main

# Directories for source and object files
SRC_DIR = ./src
OBJ_DIR = ./obj

# Source file
SRC = $(SRC_DIR)/main.f90

# Object file
OBJS = $(OBJ_DIR)/main.o

# Default rule (build the executable)
all: $(TARGET)

# Rule to link object files to create the executable, include LAPACK and BLAS
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to compile Fortran source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Clean rule to remove generated files
clean:
	rm -f $(TARGET) $(OBJ_DIR)/*.o

# Phony targets
.PHONY: all clean
