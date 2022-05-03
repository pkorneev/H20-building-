CC = gcc
CFLAGS = -std=gnu99 -Wall -Wextra -Werror -pedantic -pthread -lrt
FILES = smth.c

all:proj2

proj2: $(FILES)
	$(CC) $(CFLAGS) $(FILES) -o $@ 
