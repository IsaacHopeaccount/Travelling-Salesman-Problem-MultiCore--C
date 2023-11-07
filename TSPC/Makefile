# Variables
CC = gcc
CFLAGS = -fopenmp -O2
LDFLAGS = -fopenmp

target: cInsertion.c fInsertion.c ompcInsertion.c  coordReader.c

ci: cInsertion.c coordReader.c
	gcc cInsertion.c coordReader.c -o ci.exe -lm

fi: fInsertion.c coordReader.c
	gcc fInsertion.c coordReader.c -o ci.exe -lm