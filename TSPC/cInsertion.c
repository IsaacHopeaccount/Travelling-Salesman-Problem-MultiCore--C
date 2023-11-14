#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "coordReader.h"
#include <time.h>


void matrixDistanceCalc(double **coords, int numOfCoords, double **distanceMatrix) {
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            double x1 = coords[i][0];
            double y1 = coords[i][1];
            double x2 = coords[j][0];
            double y2 = coords[j][1];

            // Calculate the Euclidean distance
            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0;
    }
}

void cheapestInsertion(double **distanceMatrix, int numOfCoords, int *tour) {
    int visited[numOfCoords];
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
        tour[i] = -1;
    }

    // Start with the first vertex as the initial tour
    int currentVertex = 0;
    tour[0] = currentVertex;
    visited[currentVertex] = 1;

    for (int currentPosition = 1; currentPosition < numOfCoords; currentPosition++) {
        int bestVertex = -1;
        int bestPosition = -1;
        double minInsertionCost = DBL_MAX; // Initialize with a large value

        for (int vk = 0; vk < numOfCoords; vk++) {
            if (!visited[vk]) {
                for (int vn = 0; vn < currentPosition; vn++) {
                    int vi = tour[vn];
                    int viPlus1 = tour[(vn + 1) % currentPosition];
                    double insertionCost = distanceMatrix[vi][vk] + distanceMatrix[vk][viPlus1] - distanceMatrix[vi][viPlus1];

                    if (insertionCost < minInsertionCost) {
                        minInsertionCost = insertionCost;
                        bestVertex = vk;
                        bestPosition = vn;
                    }
                }
            }
        }

        for (int i = currentPosition; i > bestPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[bestPosition + 1] = bestVertex;
        visited[bestVertex] = 1;
    }
}

int main() {

    clock_t start_time = clock();

    char filename[] = "C:/Users/isaac/CLionProjects/TSPC/4096_coords.coord";

    int numOfCoords = readNumOfCoords(filename);

    if (numOfCoords < 0) {
        printf("Unable to open the file: %s\n", filename);
        return 1;
    }

    double **coords = readCoords(filename, numOfCoords);

    if (coords == NULL) {
        return 1;
    }

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    matrixDistanceCalc(coords, numOfCoords, distanceMatrix);

    int *tour = (int *)malloc(numOfCoords * sizeof(int));
    cheapestInsertion(distanceMatrix, numOfCoords, tour);

    writeTourToFile(tour,numOfCoords,"ciOut.dat");
    printf("Writing tour to data file\n");

    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", execution_time);

    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(coords);
    free(tour);

    return 0;
}
