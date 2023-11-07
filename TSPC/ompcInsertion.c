#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "coordReader.c"
#include <omp.h>

void matrixDistanceCalc(double **coords, int numOfCoords, double **distanceMatrix) {
    double start_time = omp_get_wtime(); // Capture the start time

#pragma omp parallel for collapse(2)
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            double x1 = coords[i][0];
            double y1 = coords[i][1];
            double x2 = coords[j][0];
            double y2 = coords[j][1];

            // Calculate the Euclidean distance
            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            for (int i = 0; i < numOfCoords; ++i) {
                distanceMatrix[i][i] = 0;
            }
        }
    }

    double end_time = omp_get_wtime(); // Capture the end time
    double execution_time = end_time - start_time;
    printf("Matrix calculation time: %f seconds\n", execution_time);
}

void cheapestInsertion(double **distanceMatrix, int numOfCoords, int *tour) {
    int *visited = (int *)malloc(numOfCoords * sizeof(int));
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
    }
    tour[0] = 0;
    visited[0] = 1;

    for (int currentPosition = 1; currentPosition < numOfCoords; currentPosition++) {
        int bestVertex = -1;
        int bestPosition = -1;
        double minInsertionCost = DBL_MAX; // Initialize with a large value

        #pragma omp parallel for schedule(static)
        for (int v = 0; v < numOfCoords; v++) {
            if (!visited[v]) {
                for (int i = 0; i < currentPosition; i++) {
                    int vi = tour[i];
                    int viPlus1 = tour[(i + 1) % currentPosition];
                    double insertionCost = distanceMatrix[vi][v] + distanceMatrix[v][viPlus1] - distanceMatrix[vi][viPlus1];

                    #pragma omp critical
                    {
                        if (insertionCost < minInsertionCost) {
                            minInsertionCost = insertionCost;
                            bestVertex = v;
                            bestPosition = i;
                        }
                    }
                }
            }
        }

        // Update the tour with the best insertion
        for (int i = currentPosition; i > bestPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[bestPosition + 1] = bestVertex;
        visited[bestVertex] = 1;
    }

    free(visited);
}

int main() {
    double main_start_time = omp_get_wtime(); // Capture the start time
    omp_set_num_threads(8);
    char filename[] = "C:/Users/isaac/CLionProjects/TSPMulti/16_coords.coord";

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

    printf("Final Tour: ");
    for (int i = 0; i < numOfCoords; i++) {
        printf("%d ", tour[i]);
    }
    printf("%d\n", tour[0]);

    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(coords);
    free(tour);

    double main_end_time = omp_get_wtime(); // Capture the end time
    double execution_time = main_end_time - main_start_time;
    printf("Total execution time: %f seconds\n", execution_time);

    return 0;
}
