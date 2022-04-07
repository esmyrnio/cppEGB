### Compiler & flags
CC=g++
CFLAGS=-c

### Directories
SRC_DIR=src
OBJ_DIR=build
HEAD_DIR=include

### Executable name
EXEC=EGB

### Sources (in src directory)
SOURCES=$(addprefix $(SRC_DIR)/, \
	KEH_4DEGB.cpp savgol.cpp)

### Objects (in obj directory)
OBJECTS=$(SOURCES:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

### Rules: #######################################

### General target (executable):
all: $(EXEC)

### How to make the executable:
$(EXEC): $(OBJECTS) 
	$(CC) $^ -o $@ -lm

### How to make every object:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

### How to clean up:
clean: 
	rm $(OBJECTS) $(EXEC)
