#include "includes.h"

int main(int argc, char **argv)
{
    //args
    int input_args = 6;
    if(argc < input_args){
        printf("Not enough input arguments!\n");
        return 99;
    }else if(argc > input_args){
        printf("Too many input arguments!\n");
        return 99;
    }
    int L = atoi(argv[1]);//length of the thresholds
    int min_win = atoi(argv[2]);//size of the first scale
    int max_win = atoi(argv[3]);//number of windows with the size of the last scale
    int N_iter = atoi(argv[4]);//number of simulations
    double conf_lim = atof(argv[5]);//confidence limit for thresholds
    int pol = 1;//fixing polynomial order
    //MPI inizialization
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //random vectors and dcca
    double *vec1, *vec2;
    vec1 = calloc(L, sizeof(double));
    if(!vec1){
        printf("MALLOC ERROR (vec1)\n");
        return 99;
    }
    vec2 = calloc(L, sizeof(double));
    if(!vec2){
        printf("MALLOC ERROR (vec2)\n");
        return 99;
    }
    if(rank == 0){
        printf("STARTING SIMULATION\n\n");
        system("if [ -e rho_mtx.txt ]; then rm rho_mtx.txt; fi");
        system("touch rho_mtx.txt");
    }
    //seed for random vectors
    srand48(time(0));
    for(int i = 0; i < N_iter; i++){
        if(rank == 0){
            printf("Simulation number %d / %d\r", i+1, N_iter);
            fflush(stdout);
        }
        gauss_rand_vec(L, vec1, 0, 1);
        gauss_rand_vec(L, vec2, 0, 1);
        dcca(vec1, L, vec2, L, min_win, max_win, pol, "rho.txt", rank, size);
        if(rank == 0)
            system("cat rho.txt >> rho_mtx.txt");
    }
    //free
    free(vec1); free(vec2);
    //thresholds
    if(rank == 0){
        printf("\n");
        system("rm rho.txt");
        up_down_lims("rho_mtx.txt", N_iter, conf_lim);
        system("rm rho_mtx.txt");
    }
    //finalization
    MPI_Finalize();
    return 0;
}
