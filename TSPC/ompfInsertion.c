#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "coordReader.c"
#include <omp.h>


void matrixDistanceCalc(double **coords, int numOfCoords, double **distanceMatrix) {
    double start_time = omp_get_wtime(); // Capture the start time
#pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            // Calculate the Euclidean distance
            double x1 = coords[i][0];
            double y1 = coords[i][1];
            double x2 = coords[j][0];
            double y2 = coords[j][1];
            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;

            //printf("Distance between vertex %d and vertex %d: %f\n", i, j, distance);
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0;
    }
    double end_time = omp_get_wtime(); // Capture the end time
    double execution_time = end_time - start_time;
    printf("Matrix calculation time: %f seconds\n", execution_time);
}

void farthestInsertion(double** distanceMatrix, int numOfCoords, int* tour) {
    double start_time = omp_get_wtime(); // Capture the start time

    int visited[numOfCoords];
#pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
        tour[i] = -1;
    }

    // Start with the first vertex as the initial tour
    int currentVertex = 0;
    tour[0] = currentVertex;
    visited[currentVertex] = 1;

    for (int currentPosition = 1; currentPosition < numOfCoords; currentPosition++) {
        int farthestVertex = -1;
        double maxDistance = -1.0;

        #pragma omp parallel for
        for (int vn = 0; vn < currentPosition; vn++) {
            for (int vk = 0; vk < numOfCoords; vk++) {
                if (!visited[vk]) {
                    double dist_vn_vk = distanceMatrix[tour[vn]][vk];
                    if (dist_vn_vk > maxDistance) {
                        maxDistance = dist_vn_vk;
                        farthestVertex = vk;
                    }
                }
            }
        }

        int bestPosition = -1;
        double minDistance = DBL_MAX; // Initialize with a large value

        for (int i = 0; i < currentPosition; i++) {
            int viPlus1 = tour[(i + 1) % currentPosition];
            double insertionCost = distanceMatrix[tour[i]][farthestVertex] + distanceMatrix[farthestVertex][viPlus1] - distanceMatrix[tour[i]][viPlus1];
            #pragma omp critical
            {
            if (insertionCost < minDistance) {
                minDistance = insertionCost;
                bestPosition = i;
                }
            }
        }

        for (int i = currentPosition; i > bestPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[bestPosition + 1] = farthestVertex;
        visited[farthestVertex] = 1;
    }

    double end_time = omp_get_wtime(); // Capture the end time
    double execution_time = end_time - start_time;
    printf("Insertion calculation time: %f seconds\n", execution_time);
}

int main() {
    double start_time = omp_get_wtime(); // Capture the start time

    omp_set_num_threads(12);
    char filename[] = "C:/Users/isaac/CLionProjects/TSPMulti/4096_coords.coord";

    int numOfCoords = readNumOfCoords(filename);

    if (numOfCoords < 0) {
        printf("Unable to open the file: %s\n", filename);
        return 1;
    }

    double **coords = readCoords(filename, numOfCoords);

    if (coords == NULL) {
        return 1;  // Error occurred while reading coordinates
    }

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    matrixDistanceCalc(coords, numOfCoords, distanceMatrix);

    int *tour = (int *)malloc(numOfCoords * sizeof(int));
    farthestInsertion(distanceMatrix, numOfCoords, tour);


    // Print the order of visited cities in the tour
   /* printf("Final Tour: ");
    for (int i = 0; i < numOfCoords; i++) {
        printf("%d ", tour[i]);
    }
    printf("%d\n", tour[0]);  // Print the starting city again to complete the tour
*/
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(coords);
    free(tour);

    double end_time = omp_get_wtime(); // Capture the end time
    double execution_time = end_time - start_time;
    printf("Execution time: %f seconds\n", execution_time);

    return 0;
}
