CC = gcc
CFLAGS = -Wall -Wextra 

all: Question1 Question2 Question3

Question1: Question1.c
	$(CC) $(CFLAGS) -o $@ $< -lm

Question2: Question2.c
	$(CC) $(CFLAGS) -o $@ $< -lm

Question3: Question3.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f Question1 Question2 Question3

.PHONY: all clean
