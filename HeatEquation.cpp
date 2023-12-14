#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

double time_tot = 0.1;
//double dt = 0.0002;
//double h = 0.02;
//int nPoints = 2000;
double l = 1;
double h = 0.02;
double dt = 0.0002;
double u_0 = 1;
double pi = 3.14159265358;

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int i, m;
    double time, x, sum, a;
    int nPoints = 49;

    // Calculate local domain size for each process
    //int local_nPoints = nPoints / size;
    int* local_nPoints = (int*) malloc(sizeof(int));
    *local_nPoints = nPoints / size;
    //if (rank == size - 1) {
      //  local_nPoints += nPoints % size;
    //}

    if (rank + 1 <= nPoints % size)
            (*local_nPoints) ++;

    // Initialize arrays for local data

    double* u_prev = (double*)malloc(sizeof(double) * (*local_nPoints + 2));
    double* u_next = (double*)malloc(sizeof(double) * (*local_nPoints + 2));
    //double u_prev[local_nPoints + 2];
    //double u_next[local_nPoints + 2];
    double u_exact[nPoints + 2];

    // Initialize initial conditions, including boundaries
    for (int i = 0; i <= *local_nPoints; i++) {
        u_prev[i] = u_0;
        u_next[i] = u_0;
    }
    if (rank == 0) {
        u_prev[0] = 0;
        u_next[0] = 0;
    }
    if (rank == size - 1) {
        u_prev[*local_nPoints + 1] = 0;
        u_next[*local_nPoints + 1] = 0;
    }

    // Main time loop
    //printf("For rank %d, conditions: time_tot = %f, dt = %d \n", rank, time_tot, dt);
    for (time = 0; time < time_tot; time += dt) {
        /// NEW TRANSMITION

        if (rank % 2 )
        {
            MPI_Send(&u_prev[1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        else
        {
            if (rank != size - 1)
                MPI_Recv(&u_prev[*local_nPoints + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank % 2 == 0 )
        {
            if (rank != size - 1)
                MPI_Send(&u_prev[*local_nPoints], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank % 2)
        {
            if (rank != size - 1)
                MPI_Send(&u_prev[*local_nPoints], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        else
        {
            if (rank != 0)
                MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank % 2 == 0)
        {
            if (rank != 0)
                MPI_Send(&u_prev[1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        else
        {
            if (rank != size - 1)
                MPI_Recv(&u_prev[*local_nPoints + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Implement MPI communication for exchanging boundary information
        /*if (rank > 0){
            //printf("sent to %d: %lf\n", rank-1,u_prev[1]);
            MPI_Send(&u_prev[1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (rank != size - 1){
            //printf("Before receive %d: %lf\n", rank, u_prev[local_nPoints + 1]);
            MPI_Recv(&u_prev[local_nPoints + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("received %d: %lf\n", rank, u_prev[local_nPoints + 1]);
        }
        if (rank != size - 1){
            //printf("sent to %d: %lf\n",rank+1, u_prev[local_nPoints]);
            MPI_Send(&u_prev[local_nPoints], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (rank > 0){
            MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("received %d: %lf\n",rank, u_prev[0]);
        }*/
        // Calculate the finite difference scheme for the local domain
        for (int i = 1; i <= *local_nPoints; i++) {
            u_next[i] = u_prev[i] + dt / (h * h) * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
        }
        // Update u_prev with u_next
        for (int i = 1; i <= *local_nPoints; i++)
            u_prev[i] = u_next[i];
    }
    //printf("aqui se cuelga\n");

    // Gather the results to process 0 using MPI_Gather
    // Print the values of each local array in order of their rank
    for (int process_rank = 0; process_rank < size; process_rank++) {
        if (rank == process_rank) {
            printf("Rank %d: ", rank);
            for (int i = 1; i <= *local_nPoints; i++) {
                printf("%f ", u_prev[i]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    printf("\n");

    MPI_Barrier(MPI_COMM_WORLD);
    // exact solution
    if (rank == 0) {
    printf("Exact solution: \n");
    for (i = 1; i <= nPoints; i++){
        x = i * h;
        sum = 0;
        for (m = 0; m < 5; m++){
                a =  exp(- pi * pi * (2*m+1) * (2*m+1) * time_tot) * sin( pi * (2*m+1) * x / l) / (2*m+1);
            sum = sum + 4 * u_0 * a/ pi;
        }
        u_exact[i] = sum;
        // printing the exact solution on the screen
        printf("%f  ", u_exact[i]);
    }
        printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(local_nPoints);
    free(u_prev);
    free(u_next);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
