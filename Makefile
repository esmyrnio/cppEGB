### Compiler & flags
CC=g++
CFLAGS=-std=c++17 -c -O3
LFLAGS=-lstdc++fs
### Directories
SRC_DIR=src
OBJ_DIR=build
HEAD_DIR=include

### Executable name
EXEC=EGB

### Sources (in src directory)
SOURCES=$(addprefix $(SRC_DIR)/, \
	main.cpp keh_cst.cpp eos.cpp constants.cpp savgol.cpp )

### Objects (in build directory)
OBJECTS=$(SOURCES:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

### Rules: #######################################

### General target (executable):
all: $(EXEC)

### How to make the executable:
$(EXEC): $(OBJECTS) 
	$(CC) $^ -o $@ -lm $(LFLAGS)

### How to make every object:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)

### How to clean up:
clean: 
	rm $(OBJECTS) $(EXEC)