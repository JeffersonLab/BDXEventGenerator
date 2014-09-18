CC = gcc

INC_PATH = -I./include 

SRC=$(wildcard src/*.c)
OBJ_DIR=obj
OBJS=$(patsubst src/%,$(OBJ_DIR)/%,$(SRC:.c=.o))

all: prepare $(OBJS)

$(OBJ_DIR)/%.o: src/%.c
	$(CC) -c -o $@ $^ $(INC_PATH)

clean:
	rm -f $(OBJ_DIR)/*.o

prepare:
	@mkdir -p obj

 