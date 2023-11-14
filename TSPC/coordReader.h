// coordReader.h
#ifndef COORD_READER_H
#define COORD_READER_H

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

#endif //COORD_READER_H
