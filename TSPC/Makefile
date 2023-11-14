# Variables
CC = gcc
CFLAGS = -fopenmp -O2
LDFLAGS = -fopenmp

target: cInsertion.c fInsertion.c ompcInsertion.c ompfInsertion.c coordReader.c

ci: cInsertion.c coordReader.c
	gcc cInsertion.c coordReader.h -o ci.exe -lm

fi: fInsertion.c coordReader.c
	gcc fInsertion.c coordReader.h -o fi.exe -lm

comp: ompcInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o comp.exe $^ $(LDFLAGS)  -lm

fomp: ompfInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o fomp.exe $^ $(LDFLAGS)  -lm

icomp: ompcInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o icomp.exe $^ $(LDFLAGS)  -lm

ifomp: ompfInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o ifomp.exe $^ $(LDFLAGS)  -lm
